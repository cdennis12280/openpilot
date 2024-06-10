from openpilot.common.params_pyx import Params, UnknownKeyName

class FrogPilotFunctions:
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
