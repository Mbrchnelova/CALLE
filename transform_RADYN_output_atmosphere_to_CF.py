import math
import sys


#PART 1: READ

atmosphere_name = str((sys.argv)[1])
mhdfilein_name = str((sys.argv)[2])
mhdfileout_name = str((sys.argv)[3])

f = open(atmosphere_name, "r")

lines = f.readlines()

START_READING = False
GET_NO_POINTS = False
no_radyn_points = 0
ip = -2
ir = -2
i = 1

ts_radyn = []
ds_radyn = []
hs_radyn = []

for line in lines:
	if GET_NO_POINTS:
		GET_NO_POINTS = False
		line = line.strip()
		nopoints = int(line)
		print nopoints
		ip = i	
	if line[0:6] == "* ndep":
		GET_NO_POINTS = True

	if i == ip + 2:
		START_READING = True
		ir = i
	if START_READING:
		line = line.strip()
		line = line.split("  ")
		hs_radyn.append(float(line[0]))
		ts_radyn.append(float(line[3]))
		ds_radyn.append(float(line[2]))

	if START_READING:
		if i == ir + nopoints - 1:
			START_READING = False		
	i = i + 1


f.close()





#PART 2: CONVERT

for i in range(0, nopoints):
	#cm to m
	hs_radyn[i] = hs_radyn[i]/100.0
	#g/cm3 to kg/m3
	ds_radyn[i] = ds_radyn[i]*1000.0


#PART 3: INTERPOLATE ON A NEW GRID

MHDinfile = open(mhdfilein_name ,"r")
MHDoutfile = open(mhdfileout_name, "w")


MHDlines = MHDinfile.readlines()



l = 1

lref = 1.0




#Lines below just determine the structure of the CFmesh file
getEND = False
doNOT = False

string_nodes = "!LIST_NODE"
string_data = "!LIST_STATE 1"
string_connect = "!LIST_ELEM"


idx1 = -1
idx2 = -1
idx0 = -1
i = 1

for MHDline in MHDlines:
        if MHDline[0:len(string_nodes)] == string_nodes:
                idx1 = i
        if MHDline[0:len(string_data)] == string_data:
                idx2 = i
        if MHDline[0:len(string_connect)] == string_connect:
                idx0 = i
        i = i + 1


cell_centers = []
connectivity = []
coordinates = []

nbelements = -1


#Here, insert a procedure to import your flux rope B field data
#This data could have a structure of [[x1, y1, z1, Bx1, By1, Bz1], [x2, y2, z2, Bx2, By2, Bz2], ...] 
LOUIS_DATA = []
n_LOUIS_DATA = len(LOUIS_DATA)

for MHDline in MHDlines:

	if l == 7:
                line = MHDline.split(" ")

		nbelements = int(line[1]) 
		MHDoutline = MHDline

	if getEND and MHDline[0] == "!":
		doNOT = True
		getEND = False
	if l > idx2 and MHDline[0] != "!" and MHDline[0] != "N" and len(MHDline) > 1 and doNOT == False:
		getEND = True
		vals = MHDline.split(" ")


		#We read the current cell state values
                Bx = float(vals[0])
                By = float(vals[1])
                Bz = float(vals[2])
                Ex = float(vals[3])
                Ey = float(vals[4])
                Ez = float(vals[5])
                psi = float(vals[6])
                phi = float(vals[7])
                rho0 = float(vals[8])
		rho1 = float(vals[9])
		Vx0 = float(vals[10])
		Vy0 = float(vals[11])
		Vx1 = float(vals[12])
		Vy1 = float(vals[13])
		T0 = float(vals[14])
		T1 = float(vals[15])

		cell_center = cell_centers[j]

                x = cell_center[0]
                y = cell_center[1]

		current_ionisation_fraction = rho0/rho1

		new_d = 0.0
		new_t = 0.0
		weights = 0.0

		n_mindist = -1
		mindist = 1e100
		for n in range(0, nopoints):
                        y_radyn = hs_radyn[n]  
                        d_radyn = ds_radyn[n]
                        t_radyn = ts_radyn[n] 

			distance = abs(y_radyn-y)
			if (distance == 0.0):
				distance = 1e-17

			weight = 1.0/distance
	
			new_t = new_t + t_radyn*weight
                        new_d = new_d + d_radyn*weight

			weights = weights + weight


			if distance < mindist:
				mindist = distance
				n_mindist = n

		if (weights != 0.0):
			new_t = new_t/weights
                        new_d = new_d/weights

		new_t = ts_radyn[n_mindist]
		new_d = ds_radyn[n_mindist]

		T0 = new_t
		T1 = new_t
		rho0 = new_d
		rho1 = new_d * rho0/current_ionisation_fraction

		#The final values we write into the output file
		MHDoutline = str(Bx) + " " + str(By) + " " + str(Bz) + " " + str(Ex) + " " + str(Ey) + " " + str(Ez) + " " + str(psi) + " " + str(phi) + " " + str(rho0) + " " + str(rho1) + " " + str(Vx0) + " " + str(Vy0) + " " + str(Vx1) + " " + str(Vy1) + " " + str(T0) + " " + str(T1) + "\n" 

		j = j + 1

	elif l > idx1 and l < idx2: 
		vals = MHDline.split(" ")

		#We read the node coordinaes
		x = float(vals[0])
		y = float(vals[1])

		#We manipulate the node coordinates if need be
		xMF = x * lref
		yMF = y * lref

		coordinates.append([x, y])

		#And the final node coordinates are written into the output file
		MHDoutline = str(xMF) + " " +  str(yMF) + "\n"

	elif l == idx2:
		for i in range(0, len(connectivity)):
			n1 = connectivity[i][0]
                        n2 = connectivity[i][1]
                        n3 = connectivity[i][2]
                        n4 = connectivity[i][3]

			cell_center_x = 1./4. * (coordinates[n1][0] + coordinates[n2][0] + coordinates[n3][0] + coordinates[n4][0])
                        cell_center_y = 1./4. * (coordinates[n1][1] + coordinates[n2][1] + coordinates[n3][1] + coordinates[n4][1])

			cell_centers.append([cell_center_x, cell_center_y])

		j = 0

		MHDoutline = MHDline

	elif l > idx0 and l <= idx0 + nbelements:
                vals = MHDline.split(" ")

		n1 = int(vals[0])
                n2 = int(vals[1])
                n3 = int(vals[2])
                n4 = int(vals[3])
                n5 = int(vals[4])

		connectivity.append([n1, n2, n3, n4])
		MHDoutline = MHDline
	else:
		MHDoutline = MHDline

	MHDoutfile.write(MHDoutline)

	l = l + 1

MHDinfile.close()
MHDoutfile.close()







for i in range(1, nopoints+1):

        print hs_radyn[-i]



for i in range(0, len(cell_centers)):
        if i%150 == 0:
                print cell_centers[i][1]




