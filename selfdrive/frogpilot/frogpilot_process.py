import datetime
import os
import threading

from cereal import log, messaging
from openpilot.common.params import Params
from openpilot.common.realtime import Priority, config_realtime_process
from openpilot.common.time import system_time_valid
from openpilot.system.hardware import HARDWARE

from openpilot.selfdrive.frogpilot.controls.frogpilot_planner import FrogPilotPlanner
from openpilot.selfdrive.frogpilot.controls.lib.frogpilot_functions import FrogPilotFunctions, is_url_pingable
from openpilot.selfdrive.frogpilot.controls.lib.frogpilot_variables import FrogPilotVariables
from openpilot.selfdrive.frogpilot.controls.lib.model_manager import DEFAULT_MODEL, DEFAULT_MODEL_NAME, download_model, update_models
from openpilot.selfdrive.frogpilot.controls.lib.theme_manager import ThemeManager

WIFI = log.DeviceState.NetworkType.wifi

def automatic_update_check(started, params):
  update_available = params.get_bool("UpdaterFetchAvailable")
  update_ready = params.get_bool("UpdateAvailable")
  update_state_idle = params.get("UpdaterState", encoding='utf8') == "idle"

  if update_ready and not started:
    HARDWARE.reboot()
  elif update_available:
    os.system("pkill -SIGHUP -f system.updated.updated")
  elif update_state_idle:
    os.system("pkill -SIGUSR1 -f system.updated.updated")

def time_checks(automatic_updates, deviceState, now, started, params, params_memory):
  if not is_url_pingable("https://github.com"):
    return

  screen_off = deviceState.screenBrightnessPercent == 0
  wifi_connection = deviceState.networkType == WIFI

  if screen_off and wifi_connection:
    if automatic_updates:
      automatic_update_check(started, params)

  update_maps(now, params, params_memory)
  update_models(params, params_memory, False, started)

def update_maps(now, params, params_memory):
  maps_selected = params.get("MapsSelected")

  if maps_selected is None:
    return

  day = now.day
  is_first = day == 1
  is_Sunday = now.weekday() == 6
  schedule = params.get_int("PreferredSchedule")

  maps_downloaded = os.path.exists('/data/media/0/osm/offline')

  if maps_downloaded and (schedule == 0 or (schedule == 1 and not is_Sunday) or (schedule == 2 and not is_first)):
    return

  suffix = "th" if 4 <= day <= 20 or 24 <= day <= 30 else ["st", "nd", "rd"][day % 10 - 1]
  todays_date = now.strftime(f"%B {day}{suffix}, %Y")

  if params.get("LastMapsUpdate", encoding='utf-8') == todays_date:
    return

  if params.get("OSMDownloadProgress") is None:
    params_memory.put("OSMDownloadLocations", maps_selected)
    params.put("LastMapsUpdate", todays_date)

def frogpilot_thread():
  config_realtime_process(5, Priority.CTRL_LOW)

  frogpilot_toggles = FrogPilotVariables.toggles

  params = Params()
  params_memory = Params("/dev/shm/params")

  frogpilot_functions = FrogPilotFunctions()
  frogpilot_planner = FrogPilotPlanner()
  theme_manager = ThemeManager()

  run_time_checks = False
  time_validated = system_time_valid()
  update_toggles = False

  pm = messaging.PubMaster(['frogpilotPlan'])
  sm = messaging.SubMaster(['carState', 'controlsState', 'deviceState', 'frogpilotCarControl',
                            'frogpilotCarState', 'frogpilotNavigation', 'modelV2', 'radarState'],
                            poll='modelV2', ignore_avg_freq=['radarState'])

  while True:
    sm.update()

    now = datetime.datetime.now()
    deviceState = sm['deviceState']
    started = deviceState.started

    if started and sm.updated['modelV2']:
      frogpilot_planner.update(sm['carState'], sm['controlsState'], sm['frogpilotCarControl'], sm['frogpilotCarState'],
                               sm['frogpilotNavigation'], sm['modelV2'], sm['radarState'], frogpilot_toggles)
      frogpilot_planner.publish(sm, pm, frogpilot_toggles)

    if params_memory.get("ModelToDownload", encoding='utf-8') is not None:
      download_model(params_memory)

    if FrogPilotVariables.toggles_updated:
      update_toggles = True
    elif update_toggles:
      FrogPilotVariables.update_frogpilot_params(started)

      if not frogpilot_toggles.model_manager:
        params.put("Model", DEFAULT_MODEL)
        params.put("ModelName", DEFAULT_MODEL_NAME)

      if time_validated and not started:
        frogpilot_functions.backup_toggles()

      update_toggles = False

    if now.second == 0:
      run_time_checks = True
    elif run_time_checks or not time_validated:
      time_check = threading.Thread(target=time_checks, args=(frogpilot_toggles.automatic_updates, deviceState, now, started, params, params_memory,))
      time_check.start()
      run_time_checks = False

      if not time_validated:
        time_validated = system_time_valid()
        if not time_validated:
          continue

      theme_manager.update_holiday()

def main():
  frogpilot_thread()

if __name__ == "__main__":
  main()
