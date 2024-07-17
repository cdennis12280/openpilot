import os
import requests
import shutil
import subprocess
import time
import urllib.request

from openpilot.common.basedir import BASEDIR
from openpilot.system.version import get_build_metadata

from openpilot.selfdrive.frogpilot.controls.lib.frogpilot_functions import MODELS_PATH, is_url_pingable

VERSION = "v4"

GITHUB_REPOSITORY_URL = "https://raw.githubusercontent.com/FrogAi/FrogPilot-Resources/"
GITLAB_REPOSITORY_URL = "https://gitlab.com/FrogAi/FrogPilot-Resources/-/raw/"

DEFAULT_MODEL = "north-dakota-v2"
DEFAULT_MODEL_NAME = "North Dakota V2 (Default)"

NAVIGATION_MODELS = {"certified-herbalist", "duck-amigo", "los-angeles", "recertified-herbalist"}
RADARLESS_MODELS = {"radical-turtle"}
STAGING_MODELS = {"radical-turtle", "secret-good-openpilot"}

def get_repository_url():
  if is_url_pingable("https://github.com"):
    return GITHUB_REPOSITORY_URL

  if is_url_pingable("https://gitlab.com"):
    return GITLAB_REPOSITORY_URL

  return None

def get_remote_file_size(url):
  try:
    response = requests.head(url, timeout=10)
    response.raise_for_status()
    return int(response.headers.get('Content-Length', 0))
  except requests.RequestException as request_error:
    print(f"Error fetching file size: {request_error}")
    return None

def delete_file(file):
  if os.path.exists(file):
    os.remove(file)

def handle_download_error(destination, error_message, error, params_memory):
  print(f"Error occurred: {error}")
  params_memory.put("ModelDownloadProgress", error_message)
  delete_file(destination)

def verify_download(file_path, model_url):
  if not os.path.exists(file_path):
    return False

  remote_file_size = get_remote_file_size(model_url)
  if remote_file_size is None:
    return False

  local_file_size = os.path.getsize(file_path)
  return remote_file_size == local_file_size

def download_file(destination, url, params_memory):
  try:
    with requests.get(url, stream=True, timeout=10) as r:
      r.raise_for_status()
      total_size = get_remote_file_size(url)
      downloaded_size = 0

      with open(destination, 'wb') as f:
        for chunk in r.iter_content(chunk_size=8192):
          if chunk:
            f.write(chunk)
            downloaded_size += len(chunk)
            progress = (downloaded_size / total_size) * 100
            if progress != 100:
              params_memory.put("ModelDownloadProgress", f"{progress:.0f}%")
            else:
              params_memory.put("ModelDownloadProgress", "Verifying authenticity...")

  except requests.HTTPError as http_error:
    handle_download_error(destination, f"Failed: Server error ({http_error.response.status_code})", http_error, params_memory)
  except requests.ConnectionError as connection_error:
    handle_download_error(destination, "Failed: Connection dropped...", connection_error, params_memory)
  except requests.Timeout as timeout_error:
    handle_download_error(destination, "Failed: Download timed out...", timeout_error, params_memory)
  except requests.RequestException as request_error:
    handle_download_error(destination, "Failed: Network request error. Check connection.", request_error, params_memory)
  except Exception as e:
    handle_download_error(destination, "Failed: Unexpected error.", e, params_memory)

def handle_existing_model(model, params_memory):
  error_message = f"Model {model} already exists, skipping download..."
  print(error_message)
  params_memory.put("ModelDownloadProgress", "Model already exists...")
  params_memory.remove("ModelToDownload")

def handle_verification_failure(model, model_path, model_url, params_memory):
  handle_download_error(model_path, "Issue connecting to Github, trying Gitlab", f"Model {model} verification failed. The file might be corrupted. Redownloading from Gitlab...", params_memory)
  second_repo_url = GITLAB_REPOSITORY_URL
  second_model_url = f"{second_repo_url}Models/{model}.thneed"
  download_file(model_path, second_model_url, params_memory)

  if verify_download(model_path, second_model_url):
    print(f"Model {model} redownloaded and verified successfully from Gitlab.")
  else:
    print(f"Model {model} redownload verification failed from Gitlab. The file might be corrupted.")

def download_model(params_memory):
  model = params_memory.get("ModelToDownload", encoding='utf-8')
  model_path = os.path.join(MODELS_PATH, f"{model}.thneed")

  if os.path.exists(model_path):
    handle_existing_model(model, params_memory)
    return

  repo_url = get_repository_url()
  if repo_url:
    model_url = f"{repo_url}Models/{model}.thneed"
    download_file(model_path, model_url, params_memory)

    if verify_download(model_path, model_url):
      success_message = f"Model {model} downloaded and verified successfully!"
      print(success_message)
      params_memory.put("ModelDownloadProgress", "Downloaded!")
      params_memory.remove("ModelToDownload")
    else:
      handle_verification_failure(model, model_path, model_url, params_memory)
  else:
    error_message = "Github and Gitlab are offline..."
    print(error_message)
    params_memory.put("ModelDownloadProgress", error_message)
    params_memory.remove("ModelToDownload")

def fetch_model_info(url):
  try:
    with urllib.request.urlopen(url) as response:
      model_info = [line.decode('utf-8').strip().split(' - ') for line in response.readlines()]
    return model_info
  except Exception as e:
    print(f"Failed to update models list. Error: {e}")
    return None

def update_model_lists(model_info, params):
  available_models = []
  available_model_names = []

  for model in model_info:
    model_name = model[0]
    if get_build_metadata().release_channel:
      if model_name not in STAGING_MODELS:
        available_models.append(model_name)
        available_model_names.append(model[1])
    else:
      available_models.append(model_name)
      available_model_names.append(model[1])

  params.put("AvailableModels", ','.join(available_models))
  params.put("AvailableModelsNames", ','.join(available_model_names))
  print("Models list updated successfully.")

def handle_model_deletion(params):
  model_name = params.get("ModelName", encoding='utf-8')
  if "(Default)" in model_name and model_name != DEFAULT_MODEL_NAME:
    params.put("ModelName", model_name.replace(" (Default)", ""))

  available_models = params.get("AvailableModels", encoding='utf-8').split(',')
  for model_file in os.listdir(MODELS_PATH):
    if model_file.endswith('.thneed') and model_file[:-7] not in available_models:
      delete_file(os.path.join(MODELS_PATH, model_file))
      print(f"Deleted model file: {model_file}")

def validate_current_model(params, params_memory):
  current_model = params.get("Model", encoding='utf-8')
  current_model_path = os.path.join(MODELS_PATH, f"{current_model}.thneed")
  if not os.path.exists(current_model_path):
    available_models = params.get("AvailableModels", encoding='utf-8').split(',')
    if current_model in available_models:
      params_memory.put("ModelToDownload", current_model)
    else:
      params.put("Model", DEFAULT_MODEL)
      params.put("ModelName", DEFAULT_MODEL_NAME)

def copy_default_model():
  default_model_path = os.path.join(MODELS_PATH, f"{DEFAULT_MODEL}.thneed")
  if not os.path.exists(default_model_path):
    source_path = os.path.join(BASEDIR, "selfdrive/modeld/models/supercombo.thneed")
    if os.path.exists(source_path):
      shutil.copyfile(source_path, default_model_path)
      print(f"Copied default model from {source_path} to {default_model_path}")
    else:
      print(f"Source default model not found at {source_path}. Exiting...")
      return False
  return True

def are_all_models_downloaded(params, params_memory):
  automatically_update_models = params.get_bool("AutomaticallyUpdateModels")
  available_models = params.get("AvailableModels", encoding='utf-8').split(',')
  all_models_downloaded = True

  for model in available_models:
    model_path = os.path.join(MODELS_PATH, f"{model}.thneed")
    model_url = f"{get_repository_url()}Models/{model}.thneed"

    if not os.path.exists(model_path) or not verify_download(model_path, model_url):
      if automatically_update_models:
        delete_file(model_path)
        while params_memory.get("ModelToDownload", encoding='utf-8') is not None:
          time.sleep(1)
        params_memory.put("ModelToDownload", model)
      all_models_downloaded = False

  return all_models_downloaded

def update_models(params, params_memory, boot_run=True, started=False):
  try:
    repo_url = get_repository_url()
    if repo_url is None:
      return

    url = f"{repo_url}Versions/model_names_{VERSION}.txt"
    model_info = fetch_model_info(url)
    if model_info is None:
      return

    update_model_lists(model_info, params)
    if not started and not boot_run:
      params.put_bool("ModelsDownloaded", are_all_models_downloaded(params, params_memory))
    if not boot_run:
      return

    handle_model_deletion(params)
    validate_current_model(params, params_memory)
    copy_default_model()
  except subprocess.CalledProcessError as e:
    print(f"Failed to update models. Error: {e}")
    return
