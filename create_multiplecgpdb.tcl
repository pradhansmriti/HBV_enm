mol load psf apo_HBV.psf pdb apo_HBV.pdb
mol load pdb cg_CD_avg.pdb
mol load pdb cg_AB_avg.pdb
for {set i 1} {$i < 61} {incr i} {
    puts "Iteration number: $i"
set cgCD [atomselect 1 "all"]
set cgAB [atomselect 2 "all"]
puts "C${i}D${i}"
puts "A${i}B${i}"
set selectioncd "(segid C${i} or segid D${i}) and name CA"
puts $selectioncd
set CD [atomselect 0 $selectioncd ]
$cgCD move [measure fit $cgCD $CD] 
$cgCD writepdb cg_C${i}D${i}.pdb
set selectionab "(segid A${i} or segid B${i}) and name CA"
puts $selectionab
set AB [atomselect 0 $selectionab ]
$cgAB move [measure fit $cgAB $AB] 
$cgAB writepdb cg_A${i}B${i}.pdb
mol load pdb cg_A${i}B${i}.pdb
mol load pdb cg_C${i}D${i}.pdb 
$cgCD delete
$CD delete
$cgAB delete
$AB delete
}
