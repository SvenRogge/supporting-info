import numpy as np

f_generated = "DUT-49op-DFT.xyz"
f_original = "dut49.xyz"

f_g = open(f_generated, 'r')
f_g.readline()
f_g.readline()
x_g = np.zeros(1728)
y_g = np.zeros(1728)
z_g = np.zeros(1728)
for i in xrange(1728):
    string = f_g.readline().split()
    x_g[i] = float(string[1])
    y_g[i] = float(string[2])
    z_g[i] = float(string[3])
f_g.close()
    
f_o = open(f_original, 'r')
f_o.readline()
x_o = np.zeros(1728)
y_o = np.zeros(1728)
z_o = np.zeros(1728)
for i in xrange(1728):
    string = f_o.readline().split()
    x_o[i] = float(string[2])
    y_o[i] = float(string[3])
    z_o[i] = float(string[4])
f_o.close()

line_indexes = np.zeros(1728)
number_indexes = np.zeros(1728)
for i in xrange(1728):
    # line_indexes[i] denotes for original line i, where the atom is found in the newly generated file 
    distance = (x_o[i]-x_g)**2+(y_o[i]-y_g)**2+(z_o[i]-z_g)**2
    line_indexes[i] = int(np.argmin(distance))
    # translates indexes of generated file to ordered file
    distance = (x_g[i]-x_o)**2+(y_g[i]-y_o)**2+(z_g[i]-z_o)**2
    number_indexes[i] = int(np.argmin(distance))+1


f_generated = "output.arc"
f_new = "output_ordered.arc"

f_g = open(f_generated, 'r')
f_n = open(f_new, 'w')
line = f_g.readline()
f_n.write(line)
lines=[""]*1728
for i in xrange(1728):
    lines[i]=f_g.readline()
for i in xrange(1728):
    string = lines[int(line_indexes[i])].split()
    number = number_indexes[int(string[0])-1]
    print_output = "%d %s\t%.8f\t%.8f\t%.8f\t%d" %(number, string[1], float(string[2]), float(string[3]), float(string[4]), int(string[5]))
    for item in string[6:]:
        number = number_indexes[int(item)]
        print_output += "\t%d" %number
    f_n.write(print_output+"\n")
f_n.close()

f_generated = "../output.arc"
f_new = "../dut49_cp_ordered.arc"

f_g = open(f_generated, 'r')
f_n = open(f_new, 'w')
line = f_g.readline()
f_n.write(line)
lines=[""]*1728
for i in xrange(1728):
    lines[i]=f_g.readline()
for i in xrange(1728):
    string = lines[int(line_indexes[i])].split()
    number = number_indexes[int(string[0])-1]
    print_output = "%d %s\t%.8f\t%.8f\t%.8f\t%d" %(number, string[1], float(string[2]), float(string[3]), float(string[4]), int(string[5]))
    for item in string[6:]:
        number = number_indexes[int(item)]
        print_output += "\t%d" %number
    f_n.write(print_output+"\n")
f_n.close()
