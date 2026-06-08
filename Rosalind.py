import numpy as np
import math
import time


start_time = time.perf_counter()

def parse(filename):
    str = ''
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
            Text = lines[0].strip()
            Pattern = lines[1].strip()
            return Text, Pattern
    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        return None


def PatternCount(Text, Pattern):
    count = 0
    k = len(Pattern)
    for i in range(len(Text) - k + 1):
        if Text[i:i+k] == Pattern:
            count += 1
    return count


def main(input_file):
    """
    Main function
    """
    Text, Pattern = parse(input_file)
    cnt = PatternCount(Text, Pattern)
    with open('output.txt', 'w') as output_file:
        output_file.write(f"{cnt}")
    print(f"Result saved to 'output.txt'.")


main('rosalind.txt')


end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"The function took {elapsed_time:.4f} seconds to run.")