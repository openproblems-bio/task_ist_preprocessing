import sys
from src import scs

bin_file = sys.argv[1]
image_file = sys.argv[2]

print(bin_file, flush=True)

scs.segment_cells(bin_file, image_file, align='rigid')
