import fitz  # PyMuPDF
from PIL import Image, ImageDraw, ImageEnhance, ImageFont
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
import os

# Function to convert PDF page to image with high quality
def pdf_page_to_image(pdf_path, zoom_x=2, zoom_y=2):
    pdf = fitz.open(pdf_path)
    page = pdf.load_page(0)
    matrix = fitz.Matrix(zoom_x, zoom_y)
    pix = page.get_pixmap(matrix=matrix)
    img = Image.frombytes("RGB", [pix.width, pix.height], pix.samples)
    return img

# Function to crop images
def crop_image(img, left, upper, right, lower):
    return img.crop((left, upper, right, lower))

# Function to scale both images to the same width while maintaining aspect ratio
def scale_images_to_same_width(img1, img2):
    common_width = max(img1.width, img2.width)
    img1_final_height = int((common_width / img1.width) * img1.height)
    img1_scaled = img1.resize((common_width, img1_final_height), Image.ANTIALIAS)

    img2_final_height = int((common_width / img2.width) * img2.height)
    img2_scaled = img2.resize((common_width, img2_final_height), Image.ANTIALIAS)

    return img1_scaled, img2_scaled

# Function to concatenate images vertically
def concatenate_images_vertically(img1, img2):
    width = img1.width
    total_height = img1.height + img2.height
    new_img = Image.new('RGB', (width, total_height))
    new_img.paste(img1, (0, 0))
    new_img.paste(img2, (0, img1.height))
    return new_img

# Function to draw a horizontal line
def draw_horizontal_line(image, y_position, color=(255, 0, 0), thickness=4):
    draw = ImageDraw.Draw(image)
    draw.line((0, y_position, image.width, y_position), fill=color, width=thickness)

# Function to draw a circle
def draw_circle(image, center_x, center_y, radius, color=(255, 0, 0)):
    draw = ImageDraw.Draw(image)
    left_up_point = (center_x - radius, center_y - radius)
    right_down_point = (center_x + radius, center_y + radius)
    draw.ellipse([left_up_point, right_down_point], outline=color, width=3)

# Function to add text to the image with a specified font size
def add_text(image, text, position, font_size=30, color=(255, 255, 255)):
    draw = ImageDraw.Draw(image)
    # Load a TrueType font similar to LaTeX
    font_path = "cmr10.ttf"  # Replace with the path to your Computer Modern font file
    font = ImageFont.truetype(font_path, font_size)
    # Add text to the image
    draw.text(position, text, fill=color, font=font)

# Load the images from PDFs with higher resolution (the entire page is treated as an image)
img1 = pdf_page_to_image('r22_exp.pdf', zoom_x=3, zoom_y=3)  # Adjust zoom for quality
# img2 = pdf_page_to_image('r22_num.pdf', zoom_x=3, zoom_y=3)
img2 = Image.open("r22_1.6.png").convert("RGB")


# Define cropping areas (left, upper, right, lower) for each image
crop_area_img1 = (0, 0, img1.width, img1.height / 2)  # Full height of img1
crop_area_img2 = (1050, img2.height /2 , 1410, 420)

print(img2.width/ 2, img2.height /2 , img2.width, img2.height)

# Define colors for the line and the circle
line_color = (0, 128, 0)  # Example: green color for the line
circle_color = (0, 0, 255)  # Example: blue color for the circle

# Crop the images
if img1 and img2:
    img1_cropped = crop_image(img1, *crop_area_img1)
    img2_cropped = crop_image(img2, *crop_area_img2)
    
    # Scale the cropped images to have the same width
    img1_scaled, img2_scaled = scale_images_to_same_width(img1_cropped, img2_cropped)
    
    # Combine the scaled images vertically
    combined_image = concatenate_images_vertically(img1_scaled, img2_scaled)
    
    # Draw a horizontal line between the two images
    line_y_position = img1_scaled.height  # Position for the line (at the bottom of img1)
    draw_horizontal_line(combined_image, line_y_position, color=line_color)

    # Draw a circle on the horizontal line (example parameters)
    circle_center_x = combined_image.width - 140  # Center the circle horizontally
    circle_center_y = line_y_position              # Center the circle vertically on the line
    circle_radius = 110                             # Radius of the circle
    draw_circle(combined_image, circle_center_x, circle_center_y, circle_radius, color=circle_color)

    # Add text to the combined image with a larger font size
    add_text(combined_image, "Experimental", position=(50, 10), font_size=20, color=(0, 0, 0))  # Top left in black
    add_text(combined_image, "Numerical", position=(50, 300), font_size=20, color=(0, 0, 0))  # Bottom left in black

    # Enhance the image quality
    enhancer = ImageEnhance.Contrast(combined_image)
    enhanced_image = enhancer.enhance(1.5)  # Increase contrast by a factor of 1.5

    # Save the enhanced image temporarily
    temp_image_path = "temp_enhanced_image.png"
    enhanced_image.save(temp_image_path, "PNG")

    # Save the combined image to a PDF with high quality
    pdf_filename = "r22_cropped.pdf"
    c = canvas.Canvas(pdf_filename, pagesize=letter)
    width, height = letter
    
    # Adjust the PDF size based on the image size
    img_width, img_height = enhanced_image.size
    scale_factor = min(width / img_width, height / img_height)
    new_width = int(img_width * scale_factor)
    new_height = int(img_height * scale_factor)
    
    # Draw the image onto the PDF
    c.drawImage(temp_image_path, 0, height - new_height, width=new_width, height=new_height)
    c.save()

    # Clean up the temporary image file
    os.remove(temp_image_path)

else:
    print("Error loading PDFs as images.")
