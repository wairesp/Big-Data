import lorem

# Set the desired file size in bytes
file_size = 5 * 1024**3

with open('lorem_5gb.txt', 'w') as f:
    while f.tell() < file_size:
        f.write(lorem.text())