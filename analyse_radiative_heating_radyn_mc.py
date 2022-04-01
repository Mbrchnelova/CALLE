import sys
import math

#z0, z1, d0, d1, tg1, coolt0, coolt1, heat1, advt, curt, dtnm, trl0, trl1



def writePLT(xs, ys, Qrs, Qnrs, Qtrls, Ts, rhos, mcs, filename, nnodes, coordinates, connectivity):

        f = open(filename, "w")

        line = "TITLE = Unstructured grid data\n"
        f.write(line)

        line = 'VARIABLES = "x0" "x1" "rho" "T" "Qr" "Qnr" "Qtrl" "mc" \n'
        f.write(line)

        line = 'ZONE   T= "ZONE0 Quad", N='+str(int(nnodes)) + ', E=' + str(int(len(Qrs))) + ', ZONETYPE=FEQUADRILATERAL, DATAPACKING=BLOCK, STRANDID=1, SOLUTIONTIME=0.00000000000000e+00, VARLOCATION=( [1-2]=NODAL,[3-8]=CELLCENTERED)\n'
        f.write(line)

        node_xs = []
        node_ys = []

        for i in range(0, nnodes):
                node_xs.append(coordinates[i][0])
                node_ys.append(coordinates[i][1])

        for i in range(0, len(node_xs)):
                line = str(node_xs[i])+"\n"
                f.write(line)

        for i in range(0, len(node_ys)):
                line = str(node_ys[i])+"\n"

                f.write(line)

        for i in range(0, len(rhos)-1):
		rho = (rhos[i]+rhos[i+1])/2.0
                line = str(rho)+"\n"

                f.write(line)

        for i in range(0, len(Ts)-1):
		T = (Ts[i]+Ts[i+1])/2.0
                line = str(T)+"\n"

                f.write(line)

        for i in range(0, len(Qrs)):
                line = str(Qrs[i])+"\n"

                f.write(line)

        for i in range(0, len(Qnrs)):
                line = str(Qnrs[i])+"\n"

                f.write(line)

        for i in range(0, len(Qtrls)):
                line = str(Qtrls[i])+"\n"

                f.write(line)


        for i in range(0, len(mcs)):
                line = str(mcs[i])+"\n"

                f.write(line)


        for i in range(0, len(connectivity)):
                line = str(connectivity[i][0]+1) + " " + str(connectivity[i][1]+1) + " " + str(connectivity[i][2]+1) + " " + str(connectivity[i][3]+1) + "\n"

                f.write(line)

        f.close()

        return 0









f = open("radyn_out_6e6_mc.txt", "r")

SI = False
xdim = 5e8

if SI:
	xdim = xdim * 100.


lines = f.readlines()

i = 1
ivar = -1

data = []

for i in range(0, 15):
	data.append([])


varID = -1
for line in lines:


	varstring = "Variable Data:" 

	if line[0:len(varstring)] == varstring:
		ivar = i + 1
		varID = varID + 1

	if i == ivar:
		varline = line.split()
		varline = varline[3]
		varline = varline.split(",")

		if len(varline) == 1:
			data[varID] = float(varline[0])

		else:
			varline[0] = varline[0][1:len(varline[0])]
			varline[-1] = varline[-1][0:len(varline[-1])-1]

			vararray = []
			for j in range(0, len(varline)):
				vararray.append(float(varline[j]))


			data[varID] = vararray
	i = i + 1





z1 = data[0]
d1 = data[1]
tg1 = data[2] 
coolt1 = data[3]
heat1 = data[4]
cmass1 = data[5]
trl1 = data[6]
z0 = data[7]
d0 = data[8]
coolt0 = data[9]
cmass0 = data[10]
trl0 = data[11]
dtnm = data[12] 
advt = data[13] 
curt = data[14] 


if SI:
	for i in range(0, len(z1)):
		z1[i] = z1[i]*0.01
		d1[i] = d1[i]*0.1
		tg1[i] = tg1[i]
		coolt1[i] = coolt1[i]*1e-7
		heat1[i] = heat1[i]*1e-7
		trl1[i] = trl1[i]*1e-7
		z0[i] = z0[i]*0.01
		d0[i] = d0[i]*0.1
		coolt0[i] = coolt0[i]*1e-7
		trl0[i] = trl0[i]*1e-7



theta = []

for i in range(0, len(tg1)):
	theta.append(0.55) #(5040./tg1[i])


dz0 = []
dz1 = []

for i in range(0, len(z0)-1):
	dz0.append(abs(z0[i+1] - z0[i]))

for i in range(0, len(z1)-1):
	dz1.append(abs(z1[i+1] - z1[i]))


dm0 = []
dm1 = []

for i in range(0, len(dz0)):
	dc = (d0[i] + d0[i+1])/2.0
	dm0.append(dz0[i]*dc)
for i in range(0, len(dz1)):
        dc = (d1[i] + d1[i+1])/2.0
        dm1.append(dz1[i]*dc)


dm5 = []
for i in range(0, len(dm1)):
	thetac = (theta[i] + theta[i+1])/2.0
	dm5.append(thetac*dm1[i] + (1.0-thetac)*dm0[i])


dfr1 = []
for i in range(0, len(dm5)):
	coolt1c = (coolt1[i] + coolt1[i+1])/2.0
	coolt0c = (coolt0[i] + coolt0[i+1])/2.0

	numerator = advt * coolt1c * dz1[i] + curt * coolt0c * dz0[i] 
	denominator = dm5[i] * dtnm
	dfr1.append(numerator/denominator)

dfnr1 = []
for i in range(0, len(dm5)):
	heat1c = (heat1[i] + heat1[i+1])/2.0
	dfnr1.append(-dtnm*heat1c*dz1[i]/dm5[i]/dtnm)


qr1 = []
for i in range(0, len(dfr1)):
	qr1.append(-(dfnr1[i] + dfr1[i]))


#; thin radiative losses
#qtrl1 = -( advt*trl1*dz1 + curt*trl0*dz0 )/dm5/dtnm
qtrl1 = []
for i in range(0, len(dm5)):
	trl1c = (trl1[i] + trl1[i+1])/2.0
	trl0c = (trl0[i] + trl0[i+1])/2.0
	numerator = advt * trl1c * dz1[i] + curt * trl0c * dz0[i]
	denominator = dm5[i] * dtnm
	qtrl1.append(-numerator/denominator)

for i in range(0, len(qr1)):
	print (z1[i] + z1[i+1])/2.0, (cmass1[i] + cmass1[i+1])/2.0, qr1[i], qtrl1[i], (tg1[i]+tg1[i+1])/2.0



nnodes = len(d1)  #in 1D
ncells = len(qr1)

connectivity = []
coordinates = []
xs = []
ys = []




for i in range(0, nnodes):
	coordinates.append([0.0, z1[i]])
	coordinates.append([xdim, z1[i]])
	xs.append(0.0)
	xs.append(xdim)
	ys.append(z1[i])
	ys.append(z1[i])

# 0 1
# 2 3
# 4 5
# 6 7
#

nnodes = len(coordinates) #in 2D

n = 0
for i in range(0, ncells):
	connectivity.append([n, n+1, n+3, n+2])
	n = n + 2


filename = "radyn_qr_out_CGS_6e6_mc.plt"

mcs = []
integral = 0.0

for i in range(0, ncells):
	integral = integral + dm1[i]/2.0
	mcs.append(integral)
	integral = integral + dm1[i]/2.0
	


writePLT(xs, ys, qr1, dfnr1, qtrl1, tg1, d1, mcs, filename, nnodes, coordinates, connectivity)

print len(xs), len(ys), len(qr1), len(dfnr1), len(qtrl1), len(tg1), len(d1), len(mcs), nnodes, len(coordinates), len(connectivity)

#600 600 299 299 299 300 300 299 600 600 299




