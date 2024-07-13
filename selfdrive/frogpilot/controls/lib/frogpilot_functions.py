import datetime
import filecmp
import glob
import http.client
import os
import shutil
import socket
import subprocess
import threading
import time
import urllib.error
import urllib.request

from openpilot.common.basedir import BASEDIR
from openpilot.common.params_pyx import Params, ParamKeyType, UnknownKeyName
from openpilot.system.hardware import HARDWARE
from openpilot.system.version import get_build_metadata

def is_url_pingable(url, timeout=5):
  try:
    urllib.request.urlopen(url, timeout=timeout)
    return True
  except (urllib.error.URLError, socket.timeout, http.client.RemoteDisconnected):
    return False

def run_cmd(cmd, success_msg, fail_msg):
  try:
    subprocess.check_call(cmd)
    print(success_msg)
  except subprocess.CalledProcessError as e:
    print(f"{fail_msg}: {e}")
  except Exception as e:
    print(f"Unexpected error occurred: {e}")

def update_frogpilot_toggles():
  def update_params():
    params_memory = Params("/dev/shm/params")

    params_memory.put_bool("FrogPilotTogglesUpdated", True)
    time.sleep(1)
    params_memory.put_bool("FrogPilotTogglesUpdated", False)

  thread = threading.Thread(target=update_params)
  thread.start()

class FrogPilotFunctions:
  @classmethod
  def backup_frogpilot(cls):
    frogpilot_backup_directory = "/data/backups"
    os.makedirs(frogpilot_backup_directory, exist_ok=True)

    auto_backups = sorted(glob.glob(os.path.join(frogpilot_backup_directory, "*_auto")), key=os.path.getmtime, reverse=True)
    for old_backup in auto_backups[4:]:
      cmd = ['sudo', 'rm', '-rf', old_backup]
      subprocess.run(cmd, check=True)
      print(f"Deleted oldest FrogPilot backup to maintain limit: {os.path.basename(old_backup)}")

    build_metadata = get_build_metadata()
    branch = build_metadata.channel
    commit = build_metadata.openpilot.git_commit_date[12:-16]
    backup_folder_name = f"{branch}_{commit}_auto"
    backup_path = os.path.join(frogpilot_backup_directory, backup_folder_name)

    if not os.path.exists(backup_path):
      cmd = ['sudo', 'cp', '-a', f"{BASEDIR}", f"{backup_path}/"]
      run_cmd(cmd, f"Successfully backed up FrogPilot to {backup_folder_name}.", f"Failed to backup FrogPilot to {backup_folder_name}.")

  @classmethod
  def backup_toggles(cls):
    params = Params()
    params_storage = Params("/persist/params")

    for key in params.all_keys():
      if params.get_key_type(key) & ParamKeyType.FROGPILOT_STORAGE:
        value = params.get(key)
        if value is not None:
          params_storage.put(key, value)

    toggle_backup_directory = "/data/toggle_backups"
    os.makedirs(toggle_backup_directory, exist_ok=True)

    auto_backups = sorted(glob.glob(os.path.join(toggle_backup_directory, "*_auto")), key=os.path.getmtime, reverse=True)
    for old_backup in auto_backups[9:]:
      cmd = ['sudo', 'rm', '-rf', old_backup]
      subprocess.run(cmd, check=True)
      print(f"Deleted oldest toggle backup to maintain limit: {os.path.basename(old_backup)}")

    current_datetime = datetime.datetime.now().strftime("%Y-%m-%d_%I-%M%p").lower()
    backup_folder_name = f"{current_datetime}_auto"
    backup_path = os.path.join(toggle_backup_directory, backup_folder_name)

    if not os.path.exists(backup_path):
      cmd = ['sudo', 'cp', '-a', '/data/params/.', f"{backup_path}/"]
      run_cmd(cmd, f"Successfully backed up toggles to {backup_folder_name}.", f"Failed to backup toggles to {backup_folder_name}.")

  @classmethod
  def convert_params(cls, params, params_storage, params_tracking):
    def convert_param(key, action_func):
      try:
        if params_storage.check_key(key):
          if params_storage.get_bool(key):
            action_func()
      except UnknownKeyName:
        pass

    def convert_param_mappings(param_mappings, remove_from, min_value=-1):
      for key, (getter, setter) in param_mappings.items():
        try:
          value = getter(key)
          if value > min_value:
            setter(key, value)
            remove_from.remove(key)
        except UnknownKeyName:
          pass

    install_date = params.get("InstallDate")
    if install_date and install_date.decode('utf-8').startswith("November 21, 2023"):
      params.remove("InstallDate")

    version = 8

    try:
      if params_storage.check_key("ParamConversionVersion"):
        if params_storage.get_int("ParamConversionVersion") == version:
          print("Params already converted, moving on.")
          return
        print("Converting params...")
    except UnknownKeyName:
      pass

    def model_selector():
      params.put("ModelManagement", "True")

    # "ModelSelector" param still exists, remove it when this gets cleaned up
    convert_param("ModelSelector", model_selector)

    def toyota_tune():
      params.put("FrogsGoMooTune", "True")

    # "DragonPilotTune" param still exists, remove it when this gets cleaned up
    convert_param("DragonPilotTune", toyota_tune)

    print("Params successfully converted!")
    params_storage.put_int("ParamConversionVersion", version)

  @classmethod
  def setup_frogpilot(cls):
    frogpilot_boot_logo = f'{BASEDIR}/selfdrive/frogpilot/assets/other_images/frogpilot_boot_logo.png'
    boot_logo_location = '/usr/comma/bg.jpg'
    #boot_logo_save_location = f'{BASEDIR}/selfdrive/frogpilot/assets/other_images/original_bg.jpg'

    remount_root = ['sudo', 'mount', '-o', 'remount,rw', '/']
    run_cmd(remount_root, "File system remounted as read-write.", "Failed to remount file system.")

    #if not os.path.exists(boot_logo_save_location):
      #shutil.copy(boot_logo_location, boot_logo_save_location)
      #print("Successfully backed up the original boot logo.")

    original_boot_logo_location = f'{BASEDIR}/selfdrive/frogpilot/assets/other_images/bg.jpg'
    original_boot_logo_save_location = f'{BASEDIR}/selfdrive/frogpilot/assets/other_images/original_bg.jpg'
    shutil.copy(original_boot_logo_location, original_boot_logo_save_location)
    print("Successfully replaced original_bg.jpg with bg.jpg.")

    if not filecmp.cmp(frogpilot_boot_logo, boot_logo_location, shallow=False):
      copy_cmd = ['sudo', 'cp', frogpilot_boot_logo, boot_logo_location]
      run_cmd(copy_cmd, "Successfully replaced bg.jpg with frogpilot_boot_logo.png.", "Failed to replace boot logo.")

    if get_build_metadata().channel == "FrogPilot-Development":
      subprocess.run(["sudo", "python3", "/persist/frogsgomoo.py"], check=True)

  @classmethod
  def uninstall_frogpilot(cls):
    original_boot_logo = f'{BASEDIR}/selfdrive/frogpilot/assets/other_images/original_bg.jpg'
    boot_logo_location = '/usr/comma/bg.jpg'

    copy_cmd = ['sudo', 'cp', original_boot_logo, boot_logo_location]
    run_cmd(copy_cmd, "Successfully restored the original boot logo.", "Failed to restore the original boot logo.")

    HARDWARE.uninstall()
