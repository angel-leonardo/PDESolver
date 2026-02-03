import numpy as np
from Visual import hotimage
from Read import read_vector_from_file



# --- Read Data ---
vector_data = read_vector_from_file("vector_data.bin")

if vector_data:
    print(f"Successfully read {len(vector_data)} doubles from vector_data.bin")
    image = np.array(vector_data)
    # --- Plot Image ---
    # --- Resize to 2-dim Matrix ---
    rows = int(input("Please enter in row number"))
    cols = int(input("Please enter in column number"))
    two_d_matrix = image.reshape(rows, cols)

    if (rows*cols != len(vector_data)):
        print(f"Error: Elements number mismatch!")
        input("Press ENTER to exit")
    # --- Plot ---
    hotimage(two_d_matrix)


else:
    print("Failed to read data. Exiting.")
    input("Press ENTER to exit")

