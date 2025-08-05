import json
import os
import requests
from urllib.parse import urlparse, parse_qs
from datetime import datetime

# === CONFIGURATION ===
HAR_FILE = "bhuvan.har"  
# FLOOD_DATE_CUTOFF = datetime.strptime("2019_08_01_00", "%Y_%m_%d_%H")
# BBOX_FILTER = [85.25, 20.0, 97.68, 30.55]  
OUTPUT_DIR = "flood_tiles"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def intersects_bbox(bbox_tile, bbox_filter):
    xmin, ymin, xmax, ymax = bbox_tile
    fxmin, fymin, fxmax, fymax = bbox_filter
    return not (xmax < fxmin or xmin > fxmax or ymax < fymin or ymin > fymax)

def write_world_file(wld_path, bbox, width, height):
    xmin, ymin, xmax, ymax = bbox
    pixel_x_size = (xmax - xmin) / width
    pixel_y_size = (ymin - ymax) / height  # negative for north-up
    with open(wld_path, "w") as f:
        f.write(f"{pixel_x_size}\n0.0\n0.0\n{pixel_y_size}\n{xmin}\n{ymax}\n")

# === Step 1: Load HAR file ===
with open(HAR_FILE, "r", encoding="utf-8") as f:
    har_data = json.load(f)

entries = har_data["log"]["entries"]
tile_urls = []

# === Step 2: Parse URLs ===
for entry in entries:
    url = entry["request"]["url"]
    parsed = urlparse(url)
    qs = parse_qs(parsed.query)

    # Only process WMS GetMap requests for flood layers
    if parsed.path.endswith("wms") and qs.get("REQUEST", [""])[0].lower() == "getmap":
        layer = qs.get("LAYERS", [None])[0]
        bbox_str = qs.get("BBOX", [None])[0]

        if not layer or not bbox_str or not layer.startswith("flood:"):
            continue

        try:
            # Example: flood:kl_2018_16_07 → date_str = 2018_16_07
            date_str = "_".join(layer.split(":")[1].split("_")[1:])
            flood_date = datetime.strptime(date_str, "%Y_%d_%m")  # KL format uses dd_mm
        except Exception:
            continue

        bbox = list(map(float, bbox_str.split(",")))
        tile_urls.append((url, bbox, qs, date_str))


print(f"Found {len(tile_urls)} matching flood tiles")

# === Step 3: Download ===
for i, (url, bbox, qs, date_str) in enumerate(tile_urls):
    try:
        print(f"Downloading {i+1}/{len(tile_urls)}: {url}")
        response = requests.get(url, timeout=10)
        if response.status_code != 200:
            print("Failed download.")
            continue

        image_path = os.path.join(OUTPUT_DIR, f"tile_{date_str}_{i:03d}.png")
        with open(image_path, "wb") as f:
            f.write(response.content)

        width = int(qs["WIDTH"][0])
        height = int(qs["HEIGHT"][0])
        wld_path = image_path.replace(".png", ".wld")
        write_world_file(wld_path, bbox, width, height)

    except Exception as e:
        print(f" Error: {e}")

print(f"Done. Tiles saved in: {OUTPUT_DIR}")
