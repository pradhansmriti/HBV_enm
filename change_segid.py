import shutil
for i in range(1,61):

# Source and destination file paths
    #source = "./cg_AB_avg.pdb"
    destination1 = f"cg_A{i}B{i}.pdb"
    destination2=f'cg_C{i}D{i}.pdb'
    replacements = {
    "A1": f"A{i}",
    "B1": f"B{i}"
}
    replacements2 = {
    "C1": f"C{i}",
    "D1": f"D{i}"
}
# Copy the file
    #shutil.copy(source, destination)

    #print(f"File copied from {source} to {destination}")
    #file = 'file_name.txt'
    with open(destination1, 'r') as f:
        content = f.read()
    for old, new in replacements.items():
        content = content.replace(old, new)
    
    with open(destination1, 'w') as f:
        f.write(content)
    with open(destination2, 'r') as g:
        content2 = g.read()
    for old, new in replacements2.items():
        content2 = content2.replace(old, new)
    with open(destination2, 'w') as g:
        g.write(content2)
