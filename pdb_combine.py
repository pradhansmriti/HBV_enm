import re

def get_atom_count(pdb_file):
    atom_count = 0
    with open(pdb_file, 'r') as f:
        lines_pdb=f.readlines()
        atom_count=len(lines_pdb)-2
    return atom_count

def remove_start_end_line(pdb_file):
    with open(pdb_file, 'r') as f:
        lines=f.readlines()
        return lines[1:-1]

def update_atom_number(lines_pdb,offset):
    updated_lines=[]
    print(lines_pdb)
    for linef2 in lines_pdb:
        print("updated input")
        print(linef2)
        print(linef2[6:11].strip())
        atom_num = int(linef2[6:11].strip()) + offset
        updated_lines.append(f"{linef2[:6]}{atom_num:5}{linef2[11:]}")
    return updated_lines


# Define a function to merge PDB files
def merge_pdb_files(input_file1, input_file2, output_file):
    atom_count_file1=get_atom_count(input_file1)
    print(atom_count_file1)
    with open(input_file2,'r') as f:
        line_file2=f.readlines()
        crst1=line_file2[0]
        endline=line_file2[-1]
    lines_file1=remove_start_end_line(input_file1)
    #print(lines_file1)
    lines_file2=remove_start_end_line(input_file2)
    print(lines_file2)
    updated_lines2=update_atom_number(lines_file2,atom_count_file1)
    #print(updated_lines2)
#ines = ['line1', 'line2']
#with open('filename.txt', 'w') as f:
    #f.write('\n'.join(lines))
    #lines=[crst1,lines_filw
    # Write the merged output file
    with open(output_file, 'w') as f:
        f.writelines(crst1)
        f.writelines(lines_file1)
        f.writelines(updated_lines2)
        f.write(endline)
    with open('newfile.pdb','w') as g:
           g.writelines(lines_file2)
    print(f"Merged PDB files saved as {output_file}")

# Example usage
#merge_pdb_files.get_atom_count('cg_AB_avg.pdbs')
final_file='decamer_avg.pdb'
for i in range(2,6):
    merge_pdb_files(final_file,f'cg_A{i}B{i}.pdb',final_file)
    merge_pdb_files(final_file,f'cg_C{i}D{i}.pdb',final_file)
