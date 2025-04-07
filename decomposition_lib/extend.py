import numpy as np

def extend(eY, exfactor):
    rows_eY, cols_eY = eY.shape
    # Initialize the output array with zeros
    esample = np.zeros((rows_eY * exfactor, cols_eY + exfactor - 1))
    
    # Populate the output array
    for m in range(1, exfactor + 1):
        start_row = (m - 1) * rows_eY
        end_row = m * rows_eY
        start_col = m - 1
        end_col = cols_eY + m - 1
        esample[start_row:end_row, start_col:end_col] = eY
    
    return esample
    
# Test the function if this script is run directly
if __name__ == '__main__':
    # Example matrix and expansion factor
    eY = np.array([[1, 2], [3, 4]])
    exfactor = 3
    # Call the function and print the result
    result = extend(eY, exfactor)
    print("Extended sample matrix:\n", result)
    