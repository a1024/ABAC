from skimage.io import imread
from skimage.metrics import peak_signal_noise_ratio, structural_similarity
import sys

def compare_images(image1_path, image2_path):
	# Load images
	img1 = imread(image1_path)
	img2 = imread(image2_path)

	# Check if the dimensions match
	if img1.shape != img2.shape:
		raise ValueError("Input images must have the same dimensions")

	# Compute PSNR
	psnr = peak_signal_noise_ratio(img1, img2, data_range=255)

	# Compute SSIM
	ssim = structural_similarity(img1, img2, channel_axis=-1)

	print('PSNR %12.6lf dB  SSIM %12.6lf'%(psnr, ssim))

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("Usage: python compare_images.py <image1.ppm> <image2.ppm>")
		sys.exit(1)

	compare_images(sys.argv[1], sys.argv[2])
