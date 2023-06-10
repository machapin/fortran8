import glob
import cv2
from tqdm import tqdm

# 画像から動画
file = 'v003-2/v003_e/'
video_name = 'v003-2/v003_e'
fps = 5  # [枚/s]

img_array = []
for filename in tqdm(glob.glob(file+'*.png')):
    img = cv2.imread(filename)
    height, width, layers = img.shape
    size = (width, height)
    img_array.append(img)

name = video_name+'.avi'
out = cv2.VideoWriter(name, cv2.VideoWriter_fourcc(*'DIVX'), fps, size)

for i in range(len(img_array)):
    out.write(img_array[i])
out.release()
