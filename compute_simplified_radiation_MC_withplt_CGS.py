
import math
import sys


mp = 1.6726219e-27*1000. #in CGS
kB = 1.38e-23 #in SI
pi = 3.141592654





def read_CFmesh_data(MHDinfile):

	MHDlines = MHDinfile.readlines()

	l = 1

        rhoi_arr = []
	rhon_arr = []
        Vxi_arr = []
        Vyi_arr = []
	Vxn_arr = []
	Vyn_arr = []
        Bx_arr = []
        By_arr = []
	Bz_arr = []
	Ex_arr = []
	Ey_arr = []
	Ez_arr = []
        Ti_arr = []
	Tn_arr = []
        phi_arr = []
        psi_arr = []
        h_arr = []
	w_arr = []
        mesh_arr = []


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
	cell_ws = []
	cell_hs = []
	connectivity = []
	coordinates = []

	nbelements = -1


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
		        rhoi = float(vals[8])
		        rhon = float(vals[9])
		        Vxi = float(vals[10])
		        Vyi = float(vals[11])
		        Vxn = float(vals[12])
		        Vyn = float(vals[13])
		        Bx = float(vals[0])
		        By = float(vals[1])
			Bz = float(vals[2])
			Ex = float(vals[3]) 
			Ey = float(vals[4]) 
			Ez = float(vals[5])
		        Ti = float(vals[14])
		        Tn = float(vals[15])
		        phi = float(vals[6])
		        psi = float(vals[7])

			w = cell_ws[j]
			h = cell_hs[j]
                	cell_center = cell_centers[j]

                	x = cell_center[0]
                	y = cell_center[1]

		        rhoi_arr.append(rhoi)
		        rhon_arr.append(rhon)
		        Vxi_arr.append(Vxi)
		        Vyi_arr.append(Vyi)
		        Vxn_arr.append(Vxn)
		        Vyn_arr.append(Vyn)
		        Bx_arr.append(Bx)
		        By_arr.append(By)
			Bz_arr.append(Bz)
		        Ex_arr.append(Ex)
		        Ey_arr.append(Ey)
			Ez_arr.append(Ez)
		        Ti_arr.append(Ti)
		        Tn_arr.append(Tn)
		        phi_arr.append(phi)
		        psi_arr.append(psi)


			mesh_arr.append([x, y])
			h_arr.append(h)
			w_arr.append(w)	

                	j = j + 1

        	elif l > idx1 and l < idx2:
                	vals = MHDline.split(" ")

                	#We read the node coordinaes
                	x = float(vals[0])
                	y = float(vals[1])

                	#We manipulate the node coordinates if need be
                	xMF = x 
                	yMF = y 

                	coordinates.append([x, y])


        	elif l == idx2:
        	        for i in range(0, len(connectivity)):
                	        n1 = connectivity[i][0]
                        	n2 = connectivity[i][1]
                       		n3 = connectivity[i][2]
                        	n4 = connectivity[i][3]
	
        	                cell_center_x = 1./4. * (coordinates[n1][0] + coordinates[n2][0] + coordinates[n3][0] + coordinates[n4][0])
                	        cell_center_y = 1./4. * (coordinates[n1][1] + coordinates[n2][1] + coordinates[n3][1] + coordinates[n4][1])
				cell_width = max(abs(coordinates[n4][0] - coordinates[n1][0]),abs(coordinates[n4][0] - coordinates[n2][0]), abs(coordinates[n4][0] - coordinates[n3][0]))
				cell_height = max(abs(coordinates[n4][1] - coordinates[n1][1]),abs(coordinates[n4][1] - coordinates[n2][1]), abs(coordinates[n4][1] - coordinates[n3][1]))

              	          	cell_centers.append([cell_center_x, cell_center_y])
				cell_ws.append(cell_width)
				cell_hs.append(cell_height)
                	j = 0


	        elif l > idx0 and l <= idx0 + nbelements:

        	        vals = MHDline.split(" ")

                	n1 = int(vals[0])
                	n2 = int(vals[1])
                	n3 = int(vals[2])
                	n4 = int(vals[3])
                	no = int(vals[4])
	
                	connectivity.append([n1, n2, n3, n4])
                	


        	l = l + 1

	MHDinfile.close()

	data_arr = [Bx_arr, By_arr, Bz_arr, Ex_arr, Ey_arr, Ez_arr, psi_arr, phi_arr, rhoi_arr, rhon_arr, Vxi_arr, Vyi_arr, Vxn_arr, Vyn_arr, Ti_arr, Tn_arr]

	nnodes = len(coordinates)
	
	return mesh_arr, data_arr, w_arr, h_arr, nnodes, coordinates, connectivity




def TraceRayPositive(xs, ys, ws, hs, etas, i0, j0, imax, angle, intercepted):

	stop = False
	i = i0
	j = j0
	x_travelled = 0.
	y_travelled = 0.
	it = 0
	I_mu = 0.0
	while stop == False:
		current_index = i + (j - 1)*imax - 1
		w = ws[current_index]
		h = hs[current_index]
		eta = etas[current_index]
		intercepted[current_index].append([i0, j0, angle])
		intercepted[current_index][0] = intercepted[current_index][0] + 1

		if (it == 0):
			x_travelled = 0.5*w
			eta0 = eta

		x_to_travel = w - x_travelled
		y_to_travel = h - y_travelled

		dx = x_to_travel/math.sin(angle)
		dy = y_to_travel/math.cos(angle)

		s = min(dx, dy)

		if dx < dy:
			#i crossing
			x_travelled = 0.0
			y_travelled = y_travelled + dx*math.cos(angle)
			i = i + 1

		elif dy < dx:
			#j crossing
			y_travelled = 0.0
			x_travelled = x_travelled + dy*math.sin(angle)
			j = j - 1
		else:
			#i,j crossing
			x_travelled = 0.0
			y_travelled = 0.0
			i = i + 1
			j = j - 1

		I_mu = I_mu + math.cos(angle)*s*0.5*(eta0 + eta)
 		eta0 = eta

		if (i == imax+1):
			i = i - imax
		if (j == 0):
			stop = True

		it = it + 1

	return I_mu, intercepted



def TraceRayNegative(xs, ys, ws, hs, etas, i0, j0, imax, angle, intercepted):

        stop = False
        i = i0
        j = j0
        x_travelled = 0.
        y_travelled = 0.
        it = 0
        I_mu = 0.0
        while stop == False:
                current_index = i + (j - 1)*imax - 1
                w = ws[current_index]
                h = hs[current_index]
                eta = etas[current_index]
		intercepted[current_index].append([i0, j0, -angle])
		intercepted[current_index][0] = intercepted[current_index][0] + 1
                if (it == 0):
                        x_travelled = 0.5*w
                        eta0 = eta

                x_to_travel = w - x_travelled
                y_to_travel = h - y_travelled

                dx = x_to_travel/math.sin(angle)
                dy = y_to_travel/math.cos(angle)

                s = min(dx, dy)

                if dx < dy:
                        #i crossing
                        x_travelled = 0.0
                        y_travelled = y_travelled + dx*math.cos(angle)
                        i = i - 1

                elif dy < dx:
                        #j crossing
                        y_travelled = 0.0
                        x_travelled = x_travelled + dy*math.sin(angle)
                        j = j - 1
                else:
                        #i,j crossing
                        x_travelled = 0.0
                        y_travelled = 0.0
                        i = i - 1
                        j = j - 1

                I_mu = I_mu + s*math.cos(angle)*0.5*(eta0 + eta)
                eta0 = eta

                if (i == 0):
                        i = i + imax
                if (j == 0):
                        stop = True

                it = it + 1

        return I_mu, intercepted


def TraceRayPositivePhotosphere(xs, ys, ws, hs, etas, i0, j0, imax, jmax, angle, intercepted):

        stop = False
        i = i0
        j = j0
        x_travelled = 0.
        y_travelled = 0.
        it = 0
        I_mu = 0.0
        while stop == False:
                current_index = i + (j - 1)*imax - 1
                w = ws[current_index]
                h = hs[current_index]
                eta = etas[current_index]
                intercepted[current_index].append([i0, j0, angle])
                intercepted[current_index][0] = intercepted[current_index][0] + 1

                if (it == 0):
                        x_travelled = 0.5*w
                        eta0 = eta

                x_to_travel = w - x_travelled
                y_to_travel = h - y_travelled

                dx = x_to_travel/math.sin(angle)
                dy = y_to_travel/math.cos(angle)

                s = min(dx, dy)

                if dx < dy:
                        #i crossing
                        x_travelled = 0.0
                        y_travelled = y_travelled + dx*math.cos(angle)
                        i = i + 1

                elif dy < dx:
                        #j crossing
                        y_travelled = 0.0
                        x_travelled = x_travelled + dy*math.sin(angle)
                        j = j + 1
                else:
                        #i,j crossing
                        x_travelled = 0.0
                        y_travelled = 0.0
                        i = i + 1
                        j = j + 1

                I_mu = I_mu + math.cos(angle)*s*0.5*(eta0 + eta)
                eta0 = eta

                if (i == imax+1):
                        i = i - imax
                if (j == jmax+1):
                        stop = True

                it = it + 1

        return I_mu, intercepted




def TraceRayNegativePhotosphere(xs, ys, ws, hs, etas, i0, j0, imax, jmax, angle, intercepted):

        stop = False
        i = i0
        j = j0
        x_travelled = 0.
        y_travelled = 0.
        it = 0
        I_mu = 0.0
        while stop == False:
                current_index = i + (j - 1)*imax - 1
                w = ws[current_index]
                h = hs[current_index]
                eta = etas[current_index]
                intercepted[current_index].append([i0, j0, -angle])
                intercepted[current_index][0] = intercepted[current_index][0] + 1
                if (it == 0):
                        x_travelled = 0.5*w
                        eta0 = eta

                x_to_travel = w - x_travelled
                y_to_travel = h - y_travelled

                dx = x_to_travel/math.sin(angle)
                dy = y_to_travel/math.cos(angle)

                s = min(dx, dy)

                if dx < dy:
                        #i crossing
                        x_travelled = 0.0
                        y_travelled = y_travelled + dx*math.cos(angle)
                        i = i - 1

                elif dy < dx:
                        #j crossing
                        y_travelled = 0.0
                        x_travelled = x_travelled + dy*math.sin(angle)
                        j = j + 1
                else:
                        #i,j crossing
                        x_travelled = 0.0
                        y_travelled = 0.0
                        i = i - 1
                        j = j + 1

                I_mu = I_mu + s*math.cos(angle)*0.5*(eta0 + eta)
                eta0 = eta

                if (i == 0):
                        i = i + imax
                if (j == jmax+1):
                        stop = True

                it = it + 1

        return I_mu, intercepted



def getCoronalLoss(T, ne, nH):
	Q = 0.0
	if T < 10.0**4.3:	
		Q = 0.0
	elif T < 10.0**4.6:
		Q = 10.0**(-21.85)
	elif T < 10.0**4.9:
		Q = 10.0**(-31) * T**2
	elif T < 10.0**5.4:
		Q = 10.0**(-21.2)
	elif T < 10.0**5.75:
		Q = 10.0**(-10.4) * T**(-2)
	elif T < 10.0**6.3:
		Q = 10.0**(-21.94)
	else:
		Q = 10.0**(-17.73) * T**(-2./3.)

	#print -Q*n*n
	return -Q*ne*nH



def sort_cells(imax, jmax):

	no_cells = imax * jmax

	cells_arr = []

	for j in range(1, jmax+1):
		for i in range(1, imax+1):
			
			cells_arr.append([i, j])

	return cells_arr




def computeRadiation(mesh_arr, cells_arr, data_arr, ws, hs, imax, jmax):


	Bx_arr = data_arr[0] 
	By_arr = data_arr[1]
	Bz_arr = data_arr[2]
	Ex_arr = data_arr[3]
	Ey_arr = data_arr[4]
	Ez_arr = data_arr[5]
	psi_arr = data_arr[6]
	phi_arr = data_arr[7]
	rhoi_arr = data_arr[8]
	rhon_arr = data_arr[9]
	Vxi_arr = data_arr[10]
	Vyi_arr = data_arr[11]
	Vxn_arr = data_arr[12]
	Vyn_arr = data_arr[13]
	Ti_arr = data_arr[14]
	Tn_arr = data_arr[15]

	T_arr = Ti_arr

	rho_arr = []

	for i in range(0, len(rhoi_arr)):
		rho_arr.append(rhoi_arr[i] + rhon_arr[i])


	no_cells = len(rho_arr)


	#Assumed abundances from Asplund et al. 2009
 	A_H = 10.0**(12.-12.0)
  	A_Mg = 10.0**(7.60-12.0)
  	A_Ca = 10.0**(6.34-12.0)
  	A_He = 10.0**(10.93-12.0)
	alpha = 7e-18 #1e-16

	chi_arr = []
	chi_photosphere_arr = []
	eta_arr = []
	eta_photosphere_arr = []


	print "T_c, Q_coronal, frac_HeI, chi_c"
	#A. Get chi and eta in every cell (for coronal simulations)
        for c in range(0, no_cells):
		i_c = cells_arr[c][0]
                j_c = cells_arr[c][1]
                T_c = T_arr[c]
                rho_c = rho_arr[c]
                x_c = mesh_arr[c][0]
                y_c = mesh_arr[c][1]

		#Determine the ionisation fraction of He
                g1 = 2.0*1.0*1.0
                g2 = 2.0*2.0*2.0
                dE = 1.63e-18    #J
                Saha = g2/g1 * math.exp(-(dE)/kB/T_c)
		frac_HeI = 1.- (1./(1./Saha + 1.))
		if frac_HeI > 1.0:
			frac_HeI = 1.0

		#Determine the number density of electrons, assuming that there is charge neutrality 
                rho_gcm = rho_c
		T_kK = T_c/1000.
		frac_HI = GetHIIonizationOnly(T_kK)
                neutral_fraction = frac_HI
		ni = rho_c * (1.0 - neutral_fraction) / mp
		ne = ni
		ne_cm = ne
		nH_cm = rho_c * neutral_fraction / mp


		#In actuallity, the losses are given by the local conditions in the corona, not in chromosphere
		ne_cm = 1.67e-16 / mp
		nH_cm = ne_cm
		nH_m = nH_cm * 1e6
		ne_m = ne_cm * 1e6
		#Set coronal parameters
		T_corona = 1.0e6
                Q_coronal = getCoronalLoss(T_corona, ne_m, nH_m)*10.

		alpha = 2.1e-18 

		#Determine chi --> for now for the photosphere we do the same as for the corona
                chi_c = alpha * rho_gcm * 4.407e23 * A_He * frac_HeI
		chi_photosphere = alpha * rho_gcm * 4.407e23 * A_He * frac_HeI

		#Determine eta --> for now for the photosphere we do the same as for the corona
		eta_c = -Q_coronal/4.0/pi
		eta_photosphere = -Q_coronal/4.0/pi

		if (c % 150 == 0):
			print T_c, Q_coronal, frac_HeI, chi_c 

		chi_photosphere_arr.append(chi_photosphere)
		chi_arr.append(chi_c)
		eta_photosphere_arr.append(eta_photosphere)
		eta_arr.append(eta_c)




	#B. Get tau - the optical depth at the He ionisation edge (for coronal simulations) in every cell
	tau_arr = []
	tau_photosphere_arr = []

	print "T_c, Tau at HeI"
        for c in range(0, no_cells):
                i_c = cells_arr[c][0]
                j_c = cells_arr[c][1]
                T_c = T_arr[c]
                rho_c = rho_arr[c]
                x_c = mesh_arr[c][0]
                y_c = mesh_arr[c][1]
                chi_c = chi_arr[c]

		rho_gcm = rho_c

		no_vertical_cells = jmax - j_c

		tau_c = 0.0
		tau_photosphere_c = 0.0

		y_lower = y_c
		y_upper = y_c
		chi_upper = chi_c
		chi_lower = chi_c

		#tau = int_yt^y chi(y') dy' = -int_y^yt chi(y') dy
		for v in range(1, no_vertical_cells+1):
			i_v = i_c
			j_v = j_c + v
			idx_v = i_v + (j_v-1)*imax - 1
			y_v = mesh_arr[idx_v][1]
			chi_v = chi_arr[idx_v]

			chi_upper = chi_v
			y_upper = y_v

			dy = (y_upper - y_lower)

			tau_c = tau_c + (-dy) * 0.5*(chi_upper + chi_lower)


			y_lower = y_v
			chi_lower = chi_v


		if tau_c == 0.0:
			h_c = hs[c]
			tau_c = (-h_c/2.0) * 0.5*(chi_c + chi_c)  

		tau_c = -tau_c
		tau_arr.append(tau_c)
	
		if c%150 == 0:
			print T_c, tau_c

	#B2. Get photospheric tau in every cell

	tau_photosphere_arr = []

        for c in range(0, no_cells):
                i_c = cells_arr[c][0]
                j_c = cells_arr[c][1]
                T_c = T_arr[c]
                rho_c = rho_arr[c]
                x_c = mesh_arr[c][0]
                y_c = mesh_arr[c][1]
                chi_c = chi_photosphere_arr[c]

                rho_gcm = rho_c

                no_vertical_cells = j_c

                tau_c = 0.0

		j_0 = 1
		i_0 = 1
		idx_0 = i_0 + (j_0-1)*imax - 1
                y_lower = mesh_arr[idx_0][1]
                y_upper = y_c
                chi_upper = chi_photosphere_arr[idx_0]
                chi_lower = chi_c


                #tau = int_yt^y chi(y') dy' = -int_y^yt chi(y') dy
                for v in range(1, no_vertical_cells+1):
                        i_v = i_c
                        j_v = v
                        idx_v = i_v + (j_v-1)*imax - 1
                        y_v = mesh_arr[idx_v][1]
                        chi_v = chi_photosphere_arr[idx_v]

                        chi_upper = chi_v
                        y_upper = y_v

                        dy = (y_upper - y_lower) #transform to cm

                        tau_c = tau_c + (dy) * 0.5*(chi_upper + chi_lower)

                        y_lower = y_v
                        chi_lower = chi_v

                if tau_c == 0.0:
                        h_c = hs[c]
                        tau_c = h_c/2.0 * 0.5*(chi_upper + chi_lower)

                tau_photosphere_arr.append(tau_c)




	#B3. Get total and neutral column mass in every cell
	columnmass_arr = []
	columnmass_neutral_arr = []

        for c in range(0, no_cells):
                i_c = cells_arr[c][0]
                j_c = cells_arr[c][1]
                rho_c = rho_arr[c]
                x_c = mesh_arr[c][0]
                y_c = mesh_arr[c][1]
		T_c = T_arr[c]

                no_vertical_cells = jmax - j_c

		mc_c = 0.0 
		mc_neutral_c = 0.0

		i_find = i_c
		j_find = j_c
		idx_find = i_find + (j_find-1)*imax - 1
		y_find = mesh_arr[idx_find][1]

                y_lower = y_find
                y_upper = y_find
                rho_upper = rho_arr[idx_find]
                rho_lower = rho_arr[idx_find]

		T_kK = T_c/1000.
		frac_HI = GetHIIonizationOnly(T_kK) 
                neutral_fraction = frac_HI

		rho_neutral_upper = rho_upper * neutral_fraction
		rho_neutral_lower = rho_lower * neutral_fraction

                for v in range(1, no_vertical_cells+1):
                        i_v = i_c
                        j_v = j_c + v
                        idx_v = i_v + (j_v-1)*imax - 1
                        y_v = mesh_arr[idx_v][1]
                        rho_v = rho_arr[idx_v]

                        rho_upper = rho_v
                        y_upper = y_v

                        dy = y_upper - y_lower

                        mc_c = mc_c + (dy) * 0.5*(rho_upper + rho_lower)

			T_v = T_arr[idx_v]
			T_kK = T_v/1000.
			frac_HI = GetHIIonizationOnly(T_kK) 
			neutral_fraction = frac_HI

			rho_neutral_upper = rho_upper * neutral_fraction

			mc_neutral_c = mc_neutral_c + dy * 0.5*(rho_neutral_upper + rho_neutral_lower)

                        y_lower = y_v
                        rho_lower = rho_v
			rho_neutral_lower = rho_lower * neutral_fraction

		if (mc_c == 0):
			h_c = hs[c]
			mc_c = h_c/2.0 * rho_c
			mc_neutral_c = h_c/2.0 * rho_c * neutral_fraction


                columnmass_arr.append(mc_c)
		columnmass_neutral_arr.append(mc_neutral_c)







	#C. Get local chromospheric loss in every cell
	Qrad_local_arr = []
        printidx = 0
	Qrad_H = []
	Qrad_Ca = []
	Qrad_Mg = []
	print "T_kK, logtau, math.log10(mc_gcm2), ne_cm, esc_HI, esc_MgII, esc_CaII, frac_HI, frac_MgII, frac_CaII, Q_H, Q_Mg, Q_Ca" 
	for c in range(0, no_cells):
		i_c = cells_arr[c][0]
		j_c = cells_arr[c][1]
		T_c = T_arr[c]
		rho_c = rho_arr[c]
		rhoi_c = rhoi_arr[c]
		x_c = mesh_arr[c][0]
		y_c = mesh_arr[c][1]
		T_kK = T_c/1000.


		#For magnesium and calcium, we are using the column mass to get escape probabilities, mc_gcm2 --> this is the total column mass from columnmass_arr
		mc_c = columnmass_arr[c]
		mc_gcm2 = mc_c

		#For hydrogen, we are using a depth variable determined from neutral hydrogen column mass:
		tau_c = columnmass_neutral_arr[c]

		if tau_c == 0.0:
                        tau_c = 1e-40
                logtau = math.log10(tau_c)
 
		
		Q_H, Q_Mg, Q_Ca, esc_HI, esc_MgII, esc_CaII, frac_HI, frac_MgII, frac_CaII = getFits(T_kK, logtau, mc_gcm2)


		eta = -(Q_H) / 4.0 / 3.141592654

		#In the total Q eqs, rho is the total density of the material (ions and electrons combined)
		rho_gcm = rho_c

		#The number density of the ions in the flow is the total number density times the ionisation fraction (1 - neutral fraction)
		ni = rho_gcm * (1.0 - frac_HI) / mp

		#The eletron number density is assumed to be the same as of ions for an overall charge neutrality
		ne = ni

		#We work with CGS, so no converison necessary
		ne_cm = ne


		Rad_H = -Q_H * esc_HI * frac_HI * A_H * 4.407e23 * ne_cm 
		Rad_Mg = -Q_Mg * esc_MgII * frac_MgII * A_Mg * 4.407e23 * ne_cm 
		Rad_Ca = -Q_Ca * esc_CaII * frac_CaII * A_Ca * 4.407e23 * ne_cm 

		if (c%150. == 0):
			print T_kK, logtau, math.log10(mc_gcm2), ne_cm, esc_HI, esc_MgII, esc_CaII, frac_HI, frac_MgII, frac_CaII, Q_H, Q_Mg, Q_Ca

		Qrad_local = Rad_H + Rad_Mg + Rad_Ca
		Qrad_local_arr.append(Qrad_local)

		Qrad_H.append(Rad_H)
		Qrad_Ca.append(Rad_Ca)
		Qrad_Mg.append(Rad_Mg)










	#D. Get coronal incoming radiation
	Qrad_coronal_arr = []


	#D1 --> get intensities of each ray
	ray_Is_arr = []

	intercepted = []
	for c in range(0, no_cells):
		intercepted.append([0])

	for t in range(1, imax+1):
		i_t = t
		j_t = jmax
		idx_t = i_t + (j_t-1)*imax - 1
		y_t = mesh_arr[idx_t][1]
		eta_t = eta_arr[idx_t]

		Imu_x_vert = 0.0

		y_upper = y_t
		y_lower = y_t
		eta_upper = eta_t
		eta_lower = eta_t
		I_rays = []

		for i in range(0, 9):
			dangle = pi/2.0/10.0
			angle = dangle + i*dangle
			I_ray, intercepted  = TraceRayNegative(mesh_arr[:][0], mesh_arr[:][1], ws, hs, eta_arr, i_t, j_t, imax, angle, intercepted)
			I_ray = I_ray / math.cos(angle)
			I_rays.append([-angle, I_ray])

		for i in range(0, 9):
			dangle = pi/2.0/10.0
			angle = dangle + i*dangle
			I_ray, intercepted = TraceRayPositive(mesh_arr[:][0], mesh_arr[:][1], ws, hs, eta_arr, i_t, j_t, imax, angle, intercepted)
			I_ray = I_ray / math.cos(angle)
			I_rays.append([angle, I_ray])

		ray_Is_arr.append(I_rays)

	printidx = 0

	for c in range(0, no_cells): 
		Qrad_coronal_arr.append(0.)


	for c in range(0, no_cells):
                i_c = cells_arr[c][0]
                j_c = cells_arr[c][1]
                T_c = T_arr[c]
                rho_c = rho_arr[c]
                x_c = mesh_arr[c][0]
                y_c = mesh_arr[c][1]
                tau_c = max(1e-40,tau_arr[c])
		eta_c = eta_arr[c]
		chi_c = chi_arr[c]
		
		idx_c = i_c + (j_c-1)*imax - 1

		current_rays = intercepted[c]
		
		no_rays = current_rays[0]
		cell_rays = []

		final_cell_rays = []

		angle_tol = 1e-6
		if no_rays > 0:
			for r in range(0, no_rays):
				ray = current_rays[r+1]
				ray_i0 = ray[0]
				ray_j0 = ray[1]
				ray_angle = ray[2] 
	
				for i in range(0, len(ray_Is_arr[ray_i0-1])):
					if abs(ray_angle-ray_Is_arr[ray_i0-1][i][0]) < angle_tol:
						I_mu = ray_Is_arr[ray_i0-1][i][1]
				cell_rays.append([ray_angle, I_mu])
	
			for i in range(0, len(cell_rays)):
				I_ray = cell_rays[i][1]
				angle_ray = cell_rays[i][0]
				I_tau = I_ray * math.exp(-tau_c/math.cos(angle_ray))
				cell_rays[i][1] = I_tau

			for i in range(0, len(cell_rays)):
				I_tau =  cell_rays[i][1]
				max_I = I_tau
				for j in range(0, len(cell_rays)):
					if abs(cell_rays[i][0]-cell_rays[j][0]) < angle_tol and j != i:
						if (cell_rays[j][1] > max_I):
							max_I = cell_rays[j][1]
				final_cell_rays.append([cell_rays[i][0], max_I]) 

			J2D = 0.0

			final_cell_rays = sorted(final_cell_rays, key=lambda x: x[0])
	
			for i in range(0, (len(final_cell_rays))-1):
				dangle = abs(final_cell_rays[i+1][0] - final_cell_rays[i][0])
				J2D = J2D + dangle * 0.5 * (final_cell_rays[i][1] + final_cell_rays[i+1][1])

			rho_gcm = rho_c
			Qrad_coronal = 4.0 * 3.141592654 * chi_c * J2D 
			Qrad_coronal_arr[c] = Qrad_coronal

		else:
			i_shift = 1
			not_found = True
			Q = 0.0
			while not_found:
				i_right = i_c + i_shift
				j_right = j_c
				i_left = i_c - i_shift
				j_left = j_c
				idx_right = i_right + (j_right-1)*imax - 1
				idx_left = i_left + (j_left-1)*imax - 1
				Q_right = Qrad_coronal_arr[idx_right]
				Q_left = Qrad_coronal_arr[idx_left]

				if Q_right > 0.0:
					if Q_left > 0.0:
						Q = 0.5*(Q_right + Q_left)
						not_found = False
					else:
						Q = Q_right
						not_found = False

				else:
					if Q_left > 0.0:
						Q = Q_left
						not_found = False 
				i_shift = i_shift + 1 
			Qrad_coronal_arr[c] = Q


                printidx = printidx + 1









	#E. Get photospheric incoming radiation
	Qrad_photospheric_arr = []

        for c in range(0, no_cells):
		Qrad_photospheric_arr.append(0.0)


        ray_Is_arr = []

        intercepted = []
        for c in range(0, no_cells):
                intercepted.append([0])

        for t in range(1, imax+1):
                i_t = t
                j_t = 1
                idx_t = i_t + (j_t-1)*imax - 1
                y_t = mesh_arr[idx_t][1]
                eta_t = eta_photosphere_arr[idx_t]

                Imu_x_vert = 0.0

                y_upper = y_t
                y_lower = y_t
                eta_upper = eta_t
                eta_lower = eta_t
                I_rays = []

                for i in range(0, 9):
                        dangle = pi/2.0/10.0
                        angle = dangle + i*dangle
                        I_ray, intercepted  = TraceRayNegativePhotosphere(mesh_arr[:][0], mesh_arr[:][1], ws, hs, eta_photosphere_arr, i_t, j_t, imax, jmax, angle, intercepted)
                        I_ray = I_ray / math.cos(angle)
                        I_rays.append([-angle, I_ray])

                for i in range(0, 9):
                        dangle = pi/2.0/10.0
                        angle = dangle + i*dangle
                        I_ray, intercepted = TraceRayPositivePhotosphere(mesh_arr[:][0], mesh_arr[:][1], ws, hs, eta_photosphere_arr, i_t, j_t, imax, jmax, angle, intercepted)
                        I_ray = I_ray / math.cos(angle)
                        I_rays.append([angle, I_ray])

                ray_Is_arr.append(I_rays)

        printidx = 0



        for c in range(0, no_cells):
                i_c = cells_arr[c][0]
                j_c = cells_arr[c][1]
                T_c = T_arr[c]
                rho_c = rho_arr[c]
                x_c = mesh_arr[c][0]
                y_c = mesh_arr[c][1]
                tau_p_c = max(1e-40,tau_photosphere_arr[c])
                eta_c = eta_photosphere_arr[c]
                chi_c = chi_photosphere_arr[c]

                idx_c = i_c + (j_c-1)*imax - 1

                current_rays = intercepted[c]

                no_rays = current_rays[0]
                cell_rays = []

                final_cell_rays = []

                angle_tol = 1e-6
                if no_rays > 0:
                        for r in range(0, no_rays):
                                ray = current_rays[r+1]
                                ray_i0 = ray[0]
                                ray_j0 = ray[1]
                                ray_angle = ray[2]

                                for i in range(0, len(ray_Is_arr[ray_i0-1])):
                                        if abs(ray_angle-ray_Is_arr[ray_i0-1][i][0]) < angle_tol:
                                                I_mu = ray_Is_arr[ray_i0-1][i][1]
                                cell_rays.append([ray_angle, I_mu])


                        for i in range(0, len(cell_rays)):
                                I_ray = cell_rays[i][1]
                                angle_ray = cell_rays[i][0]
                                I_tau = I_ray * math.exp(-tau_p_c/math.cos(angle_ray))
                                cell_rays[i][1] = I_tau

                        for i in range(0, len(cell_rays)):
                                I_tau =  cell_rays[i][1]
                                max_I = I_tau
                                for j in range(0, len(cell_rays)):
                                        if abs(cell_rays[i][0]-cell_rays[j][0]) < angle_tol and j != i:
                                                if (cell_rays[j][1] > max_I):
                                                        max_I = cell_rays[j][1]
                                final_cell_rays.append([cell_rays[i][0], max_I])

                        J2D = 0.0

                        final_cell_rays = sorted(final_cell_rays, key=lambda x: x[0])

                        for i in range(0, (len(final_cell_rays))-1):
                                dangle = abs(final_cell_rays[i+1][0] - final_cell_rays[i][0])
                                J2D = J2D + dangle * 0.5 * (final_cell_rays[i][1] + final_cell_rays[i+1][1])

                        rho_gcm = rho_c
                        Qrad_photospheric = 4.0 * 3.141592654 * chi_c * J2D
                        Qrad_photospheric_arr[c] = Qrad_photospheric

                else:
			print "no rays in: ", i_c, j_c
                        i_shift = 1
                        not_found = True
                        Q = 0.0
                        while not_found:
				if (i_c < imax-1):
                                	i_right = i_c + i_shift
				else:
					i_right = i_c
                                j_right = j_c

				if (i_c > 2):
                                	i_left = i_c - i_shift
				else:
					i_left = i_c
                                j_left = j_c
				#if j_c == 99:
				#	print i_c, j_c, i_shift, i_right, i_right + (j_right-1)*imax - 1
				#	print no_cells
                                idx_right = i_right + (j_right-1)*imax - 1
                                idx_left = i_left + (j_left-1)*imax - 1
                                Q_right = Qrad_photospheric_arr[idx_right]
                                Q_left = Qrad_photospheric_arr[idx_left]

                                if Q_right > 0.0:
                                        if Q_left > 0.0:
                                                Q = 0.5*(Q_right + Q_left)
                                                not_found = False
                                        else:
                                                Q = Q_right
                                                not_found = False

                                else:
                                        if Q_left > 0.0:
                                                Q = Q_left
                                                not_found = False
                                i_shift = i_shift + 1

				if i_c + i_shift == imax or i_c - i_shift == 1:
					not_found = False
					Q = 0.0
                        Qrad_photospheric_arr[c] = Q

                printidx = printidx + 1







	#F. Sum them all up
	Eloss_arr = []

	for c in range(0, no_cells):
		Qtot = Qrad_local_arr[c] + Qrad_coronal_arr[c] + Qrad_photospheric_arr[c]
		Eloss_arr.append(Qtot)

	return Eloss_arr, Qrad_local_arr, Qrad_coronal_arr, Qrad_photospheric_arr, tau_arr, columnmass_arr, Qrad_H, Qrad_Mg, Qrad_Ca









def writePLT(xs, ys, Qlocal, Qcoronal, Qphotospheric, Ts, rhos, taus, mcs,  Qrad_H, Qrad_Mg, Qrad_Ca, filename, nnodes, coordinates, connectivity):

	f = open(filename, "w")

	line = "TITLE = Unstructured grid data\n"
	f.write(line)

	line = 'VARIABLES = "x0" "x1" "rho" "T" "Qloc" "Qcor" "Qphot" "tau" "mc" "QH" "QMg" "QCa"\n'
        f.write(line)

	line = 'ZONE   T= "ZONE0 Quad", N='+str(int(nnodes)) + ', E=' + str(int(len(Ts))) + ', ZONETYPE=FEQUADRILATERAL, DATAPACKING=BLOCK, STRANDID=1, SOLUTIONTIME=0.00000000000000e+00, VARLOCATION=( [1-2]=NODAL,[3-12]=CELLCENTERED)\n'
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

        for i in range(0, len(rhos)):
                line = str(rhos[i])+"\n"

                f.write(line)

        for i in range(0, len(Ts)):
                line = str(Ts[i])+"\n"

                f.write(line)

        for i in range(0, len(Qlocal)):
                line = str(Qlocal[i])+"\n"

                f.write(line)

        for i in range(0, len(Qcoronal)):
                line = str(Qcoronal[i])+"\n"

                f.write(line)

        for i in range(0, len(Qphotospheric)):
                line = str(Qphotospheric[i])+"\n"

                f.write(line)

        for i in range(0, len(taus)):
                line = str(taus[i])+"\n"

                f.write(line)

        for i in range(0, len(mcs)):
                line = str(mcs[i])+"\n"

                f.write(line)

        for i in range(0, len(Qrad_H)):
                line = str(Qrad_H[i])+"\n"

                f.write(line)

        for i in range(0, len(Qrad_Mg)):
                line = str(Qrad_Mg[i])+"\n"

                f.write(line)

        for i in range(0, len(Qrad_Ca)):
                line = str(Qrad_Ca[i])+"\n"

                f.write(line)

	for i in range(0, len(connectivity)):
		line = str(connectivity[i][0]+1) + " " + str(connectivity[i][1]+1) + " " + str(connectivity[i][2]+1) + " " + str(connectivity[i][3]+1) + "\n"

                f.write(line)

	f.close()

	return 0


def GetHIIonizationOnly(T_kK):
  if( T_kK < 2.32844355 ) :
    frac_HI = 1.04534562

  if( 2.32844355 < T_kK and T_kK < 2.739035463 ) :
    frac_HI = 1.04534562 + (T_kK - 2.32844355) * (1.045243865 - 1.04534562) / (2.739035463 - 2.32844355)

  if( 2.739035463 < T_kK and T_kK < 3.149627377 ) :
    frac_HI = 1.045243865 + (T_kK - 2.739035463) * (1.045243865 - 1.045243865) / (3.149627377 - 2.739035463)

  if( 3.149627377 < T_kK and T_kK < 3.560219291 ) :
    frac_HI = 1.045243865 + (T_kK - 3.149627377) * (1.045243865 - 1.045243865) / (3.560219291 - 3.149627377)

  if( 3.560219291 < T_kK and T_kK < 3.970811204 ) :
    frac_HI = 1.045243865 + (T_kK - 3.560219291) * (1.045243865 - 1.045243865) / (3.970811204 - 3.560219291)

  if( 3.970811204 < T_kK and T_kK < 4.381403118 ) :
    frac_HI = 1.045243865 + (T_kK - 3.970811204) * (1.045243865 - 1.045243865) / (4.381403118 - 3.970811204)

  if( 4.381403118 < T_kK and T_kK < 4.791995032 ) :
    frac_HI = 1.045243865 + (T_kK - 4.381403118) * (1.045243865 - 1.045243865) / (4.791995032 - 4.381403118)

  if( 4.791995032 < T_kK and T_kK < 5.202586945 ) :
    frac_HI = 1.045243865 + (T_kK - 4.791995032) * (1.045243865 - 1.045243865) / (5.202586945 - 4.791995032)

  if( 5.202586945 < T_kK and T_kK < 5.613178859 ) :
    frac_HI = 1.045243865 + (T_kK - 5.202586945) * (1.045243865 - 1.045243865) / (5.613178859 - 5.202586945)

  if( 5.613178859 < T_kK and T_kK < 6.023770773 ) :
    frac_HI = 1.045243865 + (T_kK - 5.613178859) * (1.041987699 - 1.045243865) / (6.023770773 - 5.613178859)

  if( 6.023770773 < T_kK and T_kK < 6.359709611 ) :
    frac_HI = 1.041987699 + (T_kK - 6.023770773) * (1.018860202 - 1.041987699) / (6.359709611 - 6.023770773)

  if( 6.359709611 < T_kK and T_kK < 6.565005568 ) :
    frac_HI = 1.018860202 + (T_kK - 6.359709611) * (0.9923566131 - 1.018860202) / (6.565005568 - 6.359709611)

  if( 6.565005568 < T_kK and T_kK < 6.714311718 ) :
    frac_HI = 0.9923566131 + (T_kK - 6.565005568) * (0.9638142866 - 0.9923566131) / (6.714311718 - 6.565005568)

  if( 6.714311718 < T_kK and T_kK < 6.863617869 ) :
    frac_HI = 0.9638142866 + (T_kK - 6.714311718) * (0.9361114403 - 0.9638142866) / (6.863617869 - 6.714311718)

  if( 6.863617869 < T_kK and T_kK < 6.99426075 ) :
    frac_HI = 0.9361114403 + (T_kK - 6.863617869) * (0.9094346254 - 0.9361114403) / (6.99426075 - 6.863617869)

  if( 6.99426075 < T_kK and T_kK < 7.106240363 ) :
    frac_HI = 0.9094346254 + (T_kK - 6.99426075) * (0.8855560777 - 0.9094346254) / (7.106240363 - 6.99426075)

  if( 7.106240363 < T_kK and T_kK < 7.218219976 ) :
    frac_HI = 0.8855560777 + (T_kK - 7.106240363) * (0.8616775301 - 0.8855560777) / (7.218219976 - 7.106240363)

  if( 7.218219976 < T_kK and T_kK < 7.348862858 ) :
    frac_HI = 0.8616775301 + (T_kK - 7.218219976) * (0.8362132976 - 0.8616775301) / (7.348862858 - 7.218219976)

  if( 7.348862858 < T_kK and T_kK < 7.516832277 ) :
    frac_HI = 0.8362132976 + (T_kK - 7.348862858) * (0.8133234711 - 0.8362132976) / (7.516832277 - 7.348862858)
  if( 7.516832277 < T_kK and T_kK < 7.722128234 ) :
    frac_HI = 0.8133234711 + (T_kK - 7.516832277) * (0.7864974153 - 0.8133234711) / (7.722128234 - 7.516832277)

  if( 7.722128234 < T_kK and T_kK < 8.039403803 ) :
    frac_HI = 0.7864974153 + (T_kK - 7.722128234) * (0.7620592142 - 0.7864974153) / (8.039403803 - 7.722128234)

  if( 8.039403803 < T_kK and T_kK < 8.449995717 ) :
    frac_HI = 0.7620592142 + (T_kK - 8.039403803) * (0.7415046689 - 0.7620592142) / (8.449995717 - 8.039403803)

  if( 8.449995717 < T_kK and T_kK < 8.860587631 ) :
    frac_HI = 0.7415046689 + (T_kK - 8.449995717) * (0.7263431479 - 0.7415046689) / (8.860587631 - 8.449995717)

  if( 8.860587631 < T_kK and T_kK < 9.271179544 ) :
    frac_HI = 0.7263431479 + (T_kK - 8.860587631) * (0.7094517889 - 0.7263431479) / (9.271179544 - 8.860587631)

  if( 9.271179544 < T_kK and T_kK < 9.681771458 ) :
    frac_HI = 0.7094517889 + (T_kK - 9.271179544) * (0.6929674506 - 0.7094517889) / (9.681771458 - 9.271179544)

  if( 9.681771458 < T_kK and T_kK < 10.09236337 ) :
    frac_HI = 0.6929674506 + (T_kK - 9.681771458) * (0.6763813571 - 0.6929674506) / (10.09236337 - 9.681771458)

  if( 10.09236337 < T_kK and T_kK < 10.50295529 ) :
    frac_HI = 0.6763813571 + (T_kK - 10.09236337) * (0.6592864878 - 0.6763813571) / (10.50295529 - 10.09236337)

  if( 10.50295529 < T_kK and T_kK < 10.9135472 ) :
    frac_HI = 0.6592864878 + (T_kK - 10.50295529) * (0.6429039047 - 0.6592864878) / (10.9135472 - 10.50295529)

  if( 10.9135472 < T_kK and T_kK < 11.32413911 ) :
    frac_HI = 0.6429039047 + (T_kK - 10.9135472) * (0.6267248319 - 0.6429039047) / (11.32413911 - 10.9135472)

  if( 11.32413911 < T_kK and T_kK < 11.73473103 ) :
    frac_HI = 0.6267248319 + (T_kK - 11.32413911) * (0.6094264523 - 0.6267248319) / (11.73473103 - 11.32413911)

  if( 11.73473103 < T_kK and T_kK < 12.14532294 ) :
    frac_HI = 0.6094264523 + (T_kK - 11.73473103) * (0.5917210519 - 0.6094264523) / (12.14532294 - 11.73473103)

  if( 12.14532294 < T_kK and T_kK < 12.55591485 ) :
    frac_HI = 0.5917210519 + (T_kK - 12.14532294) * (0.5732016101 - 0.5917210519) / (12.55591485 - 12.14532294)

  if( 12.55591485 < T_kK and T_kK < 12.96650677 ) :
    frac_HI = 0.5732016101 + (T_kK - 12.55591485) * (0.5526470648 - 0.5732016101) / (12.96650677 - 12.55591485)

  if( 12.96650677 < T_kK and T_kK < 13.37709868 ) :
    frac_HI = 0.5526470648 + (T_kK - 12.96650677) * (0.5286328436 - 0.5526470648) / (13.37709868 - 12.96650677)

  if( 13.37709868 < T_kK and T_kK < 13.78769059 ) :
    frac_HI = 0.5286328436 + (T_kK - 13.37709868) * (0.5050256431 - 0.5286328436) / (13.78769059 - 13.37709868)

  if( 13.78769059 < T_kK and T_kK < 14.19828251 ) :
    frac_HI = 0.5050256431 + (T_kK - 13.78769059) * (0.4811131771 - 0.5050256431) / (14.19828251 - 13.78769059)

  if( 14.19828251 < T_kK and T_kK < 14.60887442 ) :
    frac_HI = 0.4811131771 + (T_kK - 14.19828251) * (0.4575059766 - 0.4811131771) / (14.60887442 - 14.19828251)
  if( 14.60887442 < T_kK and T_kK < 15.01946634 ) :
    frac_HI = 0.4575059766 + (T_kK - 14.60887442) * (0.4340005313 - 0.4575059766) / (15.01946634 - 14.60887442)

  if( 15.01946634 < T_kK and T_kK < 15.43005825 ) :
    frac_HI = 0.4340005313 + (T_kK - 15.01946634) * (0.4103933307 - 0.4340005313) / (15.43005825 - 15.01946634)

  if( 15.43005825 < T_kK and T_kK < 15.84065016 ) :
    frac_HI = 0.4103933307 + (T_kK - 15.43005825) * (0.3903475613 - 0.4103933307) / (15.84065016 - 15.43005825)

  if( 15.84065016 < T_kK and T_kK < 16.25124208 ) :
    frac_HI = 0.3903475613 + (T_kK - 15.84065016) * (0.3705053023 - 0.3903475613) / (16.25124208 - 15.84065016)

  if( 16.25124208 < T_kK and T_kK < 16.64317072 ) :
    frac_HI = 0.3705053023 + (T_kK - 16.25124208) * (0.3467963466 - 0.3705053023) / (16.64317072 - 16.25124208)

  if( 16.64317072 < T_kK and T_kK < 17.0164361 ) :
    frac_HI = 0.3467963466 + (T_kK - 16.64317072) * (0.3218358023 - 0.3467963466) / (17.0164361 - 16.64317072)

  if( 17.0164361 < T_kK and T_kK < 17.3710382 ) :
    frac_HI = 0.3218358023 + (T_kK - 17.0164361) * (0.2980443119 - 0.3218358023) / (17.3710382 - 17.0164361)

  if( 17.3710382 < T_kK and T_kK < 17.70697704 ) :
    frac_HI = 0.2980443119 + (T_kK - 17.3710382) * (0.2730464573 - 0.2980443119) / (17.70697704 - 17.3710382)

  if( 17.70697704 < T_kK and T_kK < 18.02425261 ) :
    frac_HI = 0.2730464573 + (T_kK - 17.70697704) * (0.2488569911 - 0.2730464573) / (18.02425261 - 17.70697704)

  if( 18.02425261 < T_kK and T_kK < 18.32286491 ) :
    frac_HI = 0.2488569911 + (T_kK - 18.02425261) * (0.2259111992 - 0.2488569911) / (18.32286491 - 18.02425261)

  if( 18.32286491 < T_kK and T_kK < 18.65880375 ) :
    frac_HI = 0.2259111992 + (T_kK - 18.32286491) * (0.201622239 - 0.2259111992) / (18.65880375 - 18.32286491)

  if( 18.65880375 < T_kK and T_kK < 19.03206913 ) :
    frac_HI = 0.201622239 + (T_kK - 18.65880375) * (0.1783406551 - 0.201622239) / (19.03206913 - 18.65880375)

  if( 19.03206913 < T_kK and T_kK < 19.42399777 ) :
    frac_HI = 0.1783406551 + (T_kK - 19.03206913) * (0.1542755563 - 0.1783406551) / (19.42399777 - 19.03206913)

  if( 19.42399777 < T_kK and T_kK < 19.83458969 ) :
    frac_HI = 0.1542755563 + (T_kK - 19.42399777) * (0.131380642 - 0.1542755563) / (19.83458969 - 19.42399777)

  if( 19.83458969 < T_kK and T_kK < 20.2451816 ) :
    frac_HI = 0.131380642 + (T_kK - 19.83458969) * (0.1158121003 - 0.131380642) / (20.2451816 - 19.83458969)

  if( 20.2451816 < T_kK and T_kK < 20.65577351 ) :
    frac_HI = 0.1158121003 + (T_kK - 20.2451816) * (0.1006505793 - 0.1158121003) / (20.65577351 - 20.2451816)

  if( 20.65577351 < T_kK and T_kK < 21.06636543 ) :
    frac_HI = 0.1006505793 + (T_kK - 20.65577351) * (0.08660836517 - 0.1006505793) / (21.06636543 - 20.65577351)

  if( 21.06636543 < T_kK and T_kK < 21.47695734 ) :
    frac_HI = 0.08660836517 + (T_kK - 21.06636543) * (0.07643284771 - 0.08660836517) / (21.47695734 - 21.06636543)

  if( 21.47695734 < T_kK and T_kK < 21.88754925 ) :
    frac_HI = 0.07643284771 + (T_kK - 21.47695734) * (0.0668678613 - 0.07643284771) / (21.88754925 - 21.47695734)
  if( 21.88754925 < T_kK and T_kK < 22.29814117 ) :
    frac_HI = 0.0668678613 + (T_kK - 21.88754925) * (0.05730287489 - 0.0668678613) / (22.29814117 - 21.88754925)

  if( 22.29814117 < T_kK and T_kK < 22.70873308 ) :
    frac_HI = 0.05730287489 + (T_kK - 22.29814117) * (0.04987474714 - 0.05730287489) / (22.70873308 - 22.29814117)

  if( 22.70873308 < T_kK and T_kK < 23.119325 ) :
    frac_HI = 0.04987474714 + (T_kK - 22.70873308) * (0.04315890562 - 0.04987474714) / (23.119325 - 22.70873308)

  if( 23.119325 < T_kK and T_kK < 23.52991691 ) :
    frac_HI = 0.04315890562 + (T_kK - 23.119325) * (0.0368500848 - 0.04315890562) / (23.52991691 - 23.119325)

  if( 23.52991691 < T_kK and T_kK < 23.94050882 ) :
    frac_HI = 0.0368500848 + (T_kK - 23.52991691) * (0.03105003984 - 0.0368500848) / (23.94050882 - 23.52991691)

  if( 23.94050882 < T_kK and T_kK < 24.35110074 ) :
    frac_HI = 0.03105003984 + (T_kK - 23.94050882) * (0.02708158803 - 0.03105003984) / (24.35110074 - 23.94050882)

  if( 24.35110074 < T_kK and T_kK < 24.76169265 ) :
    frac_HI = 0.02708158803 + (T_kK - 24.35110074) * (0.02311313623 - 0.02708158803) / (24.76169265 - 24.35110074)

  if( 24.76169265 < T_kK and T_kK < 25.17228456 ) :
    frac_HI = 0.02311313623 + (T_kK - 24.76169265) * (0.01924643959 - 0.02311313623) / (25.17228456 - 24.76169265)

  if( 25.17228456 < T_kK and T_kK < 25.58287648 ) :
    frac_HI = 0.01924643959 + (T_kK - 25.17228456) * (0.01629553953 - 0.01924643959) / (25.58287648 - 25.17228456)

  if( 25.58287648 < T_kK and T_kK < 25.99346839 ) :
    frac_HI = 0.01629553953 + (T_kK - 25.58287648) * (0.01497272226 - 0.01629553953) / (25.99346839 - 25.58287648)

  if( 25.99346839 < T_kK and T_kK < 26.4040603 ) :
    frac_HI = 0.01497272226 + (T_kK - 25.99346839) * (0.01324288429 - 0.01497272226) / (26.4040603 - 25.99346839)

  if( 26.4040603 < T_kK and T_kK < 26.81465222 ) :
    frac_HI = 0.01324288429 + (T_kK - 26.4040603) * (0.01151304632 - 0.01324288429) / (26.81465222 - 26.4040603)

  if( 26.81465222 < T_kK and T_kK < 27.22524413 ) :
    frac_HI = 0.01151304632 + (T_kK - 26.81465222) * (0.01019022905 - 0.01151304632) / (27.22524413 - 26.81465222)

  if( 27.22524413 < T_kK and T_kK < 27.63583605 ) :
    frac_HI = 0.01019022905 + (T_kK - 27.22524413) * (0.009579698004 - 0.01019022905) / (27.63583605 - 27.22524413)

  if( 27.63583605 < T_kK and T_kK < 28.04642796 ) :
    frac_HI = 0.009579698004 + (T_kK - 27.63583605) * (0.008765656607 - 0.009579698004) / (28.04642796 - 27.63583605)

  if( 28.04642796 < T_kK and T_kK < 28.45701987 ) :
    frac_HI = 0.008765656607 + (T_kK - 28.04642796) * (0.007748104861 - 0.008765656607) / (28.45701987 - 28.04642796)

  if( 28.45701987 < T_kK and T_kK < 28.86761179 ) :
    frac_HI = 0.007748104861 + (T_kK - 28.45701987) * (0.006527042766 - 0.007748104861) / (28.86761179 - 28.45701987)

  if( 28.86761179 < T_kK and T_kK < 29.2782037 ) :
    frac_HI = 0.006527042766 + (T_kK - 28.86761179) * (0.006527042766 - 0.006527042766) / (29.2782037 - 28.86761179)

  if( 29.2782037 < T_kK and T_kK < 29.68879561 ) :
    frac_HI = 0.006527042766 + (T_kK - 29.2782037) * (0.005407735845 - 0.006527042766) / (29.68879561 - 29.2782037)

  if( 29.68879561 < T_kK and T_kK < 29.95008138 ) :
    frac_HI = 0.005407735845 + (T_kK - 29.68879561) * (0.005407735845 - 0.005407735845) / (29.95008138 - 29.68879561)

  if( T_kK > 29.95008138 ) :
    frac_HI = 0.0; #0.005407735845

  if frac_HI > 1.0:
	frac_HI =1.0

  return frac_HI








def getFits(T_kK, logtau, mc):


  if( T_kK < 7.692503177 ) : 
    log10Q_HI_Lyman = -24.98152317 
    
  if( 7.692503177 < T_kK and T_kK < 7.977128335 ) : 
    log10Q_HI_Lyman = -24.98152317 + (T_kK - 7.692503177) * (-24.72033061 - -24.98152317) / (7.977128335 - 7.692503177) 
    
  if( 7.977128335 < T_kK and T_kK < 8.297331639 ) : 
    log10Q_HI_Lyman = -24.72033061 + (T_kK - 7.977128335) * (-24.45911855 - -24.72033061) / (8.297331639 - 7.977128335) 
    
  if( 8.297331639 < T_kK and T_kK < 8.653113088 ) : 
    log10Q_HI_Lyman = -24.45911855 + (T_kK - 8.297331639) * (-24.19788697 - -24.45911855) / (8.653113088 - 8.297331639) 
    
  if( 8.653113088 < T_kK and T_kK < 9.115628971 ) : 
    log10Q_HI_Lyman = -24.19788697 + (T_kK - 8.653113088) * (-23.92124177 - -24.19788697) / (9.115628971 - 8.653113088) 
    
  if( 9.115628971 < T_kK and T_kK < 9.684879288 ) : 
    log10Q_HI_Lyman = -23.92124177 + (T_kK - 9.115628971) * (-23.61382787 - -23.92124177) / (9.684879288 - 9.115628971) 
    
  if( 9.684879288 < T_kK and T_kK < 10.28970775 ) : 
    log10Q_HI_Lyman = -23.61382787 + (T_kK - 9.684879288) * (-23.29103937 - -23.61382787) / (10.28970775 - 9.684879288) 
    
  if( 10.28970775 < T_kK and T_kK < 10.9656925 ) : 
    log10Q_HI_Lyman = -23.29103937 + (T_kK - 10.28970775) * (-22.99892202 - -23.29103937) / (10.9656925 - 10.28970775) 
    
  if( 10.9656925 < T_kK and T_kK < 11.71283355 ) : 
    log10Q_HI_Lyman = -22.99892202 + (T_kK - 10.9656925) * (-22.69141057 - -22.99892202) / (11.71283355 - 10.9656925) 
    
  if( 11.71283355 < T_kK and T_kK < 12.49555273 ) : 
    log10Q_HI_Lyman = -22.69141057 + (T_kK - 11.71283355) * (-22.44529994 - -22.69141057) / (12.49555273 - 11.71283355) 
    
  if( 12.49555273 < T_kK and T_kK < 13.52731893 ) : 
    log10Q_HI_Lyman = -22.44529994 + (T_kK - 12.49555273) * (-22.12227731 - -22.44529994) / (13.52731893 - 12.49555273) 
    
  if( 13.52731893 < T_kK and T_kK < 14.66581957 ) : 
    log10Q_HI_Lyman = -22.12227731 + (T_kK - 13.52731893) * (-21.84526141 - -22.12227731) / (14.66581957 - 13.52731893) 
    
  if( 14.66581957 < T_kK and T_kK < 15.62642948 ) : 
    log10Q_HI_Lyman = -21.84526141 + (T_kK - 14.66581957) * (-21.61440832 - -21.84526141) / (15.62642948 - 14.66581957) 
    
  if( 15.62642948 < T_kK and T_kK < 16.94282084 ) : 
    log10Q_HI_Lyman = -21.61440832 + (T_kK - 15.62642948) * (-21.36800503 - -21.61440832) / (16.94282084 - 15.62642948) 
    
  if( 16.94282084 < T_kK and T_kK < 18.40152478 ) : 
    log10Q_HI_Lyman = -21.36800503 + (T_kK - 16.94282084) * (-21.13687879 - -21.36800503) / (18.40152478 - 16.94282084) 
    
  if( 18.40152478 < T_kK and T_kK < 19.78907243 ) : 
    log10Q_HI_Lyman = -21.13687879 + (T_kK - 18.40152478) * (-20.93650174 - -21.13687879) / (19.78907243 - 18.40152478) 
    
  if( 19.78907243 < T_kK and T_kK < 22.20838628 ) : 
    log10Q_HI_Lyman = -20.93650174 + (T_kK - 19.78907243) * (-20.62807327 - -20.93650174) / (22.20838628 - 19.78907243) 
    
  if( 22.20838628 < T_kK and T_kK < 23.88055909 ) : 
    log10Q_HI_Lyman = -20.62807327 + (T_kK - 22.20838628) * (-20.4736054 - -20.62807327) / (23.88055909 - 22.20838628) 
    
  if( 23.88055909 < T_kK and T_kK < 25.73062262 ) : 
    log10Q_HI_Lyman = -20.4736054 + (T_kK - 23.88055909) * (-20.30368488 - -20.4736054) / (25.73062262 - 23.88055909) 
    
  if( 25.73062262 < T_kK and T_kK < 27.50952986 ) : 
    log10Q_HI_Lyman = -20.30368488 + (T_kK - 25.73062262) * (-20.16451356 - -20.30368488) / (27.50952986 - 25.73062262) 
    
  if( 27.50952986 < T_kK and T_kK < 29.14612452 ) : 
    log10Q_HI_Lyman = -20.16451356 + (T_kK - 27.50952986) * (-20.04077536 - -20.16451356) / (29.14612452 - 27.50952986) 
    
  if( 29.14612452 < T_kK and T_kK < 29.96442186 ) : 
    log10Q_HI_Lyman = -20.04077536 + (T_kK - 29.14612452) * (-19.97890627 - -20.04077536) / (29.96442186 - 29.14612452) 
    
  if( T_kK > 29.96442186 ) : 
    log10Q_HI_Lyman = -19.97890627 
    
  if( T_kK < 4.723021878 ) : 
    log10Q_HI_neutral = -25.00763359 
    
  if( 4.723021878 < T_kK and T_kK < 5.082934444 ) : 
    log10Q_HI_neutral = -25.00763359 + (T_kK - 4.723021878) * (-24.77862595 - -25.00763359) / (5.082934444 - 4.723021878) 
    
  if( 5.082934444 < T_kK and T_kK < 5.380707096 ) : 
    log10Q_HI_neutral = -24.77862595 + (T_kK - 5.082934444) * (-24.73282443 - -24.77862595) / (5.380707096 - 5.082934444) 
    
  if( 5.380707096 < T_kK and T_kK < 5.798790786 ) : 
    log10Q_HI_neutral = -24.73282443 + (T_kK - 5.380707096) * (-24.58778626 - -24.73282443) / (5.798790786 - 5.380707096) 
    
  if( 5.798790786 < T_kK and T_kK < 6.335597997 ) : 
    log10Q_HI_neutral = -24.58778626 + (T_kK - 5.798790786) * (-24.45038168 - -24.58778626) / (6.335597997 - 5.798790786) 
    
  if( 6.335597997 < T_kK and T_kK < 6.752434352 ) : 
    log10Q_HI_neutral = -24.45038168 + (T_kK - 6.335597997) * (-24.38931298 - -24.45038168) / (6.752434352 - 6.335597997) 
    
  if( 6.752434352 < T_kK and T_kK < 7.229142741 ) : 
    log10Q_HI_neutral = -24.38931298 + (T_kK - 6.752434352) * (-24.29770992 - -24.38931298) / (7.229142741 - 6.752434352) 
    
  if( 7.229142741 < T_kK and T_kK < 7.587014215 ) : 
    log10Q_HI_neutral = -24.29770992 + (T_kK - 7.229142741) * (-24.20610687 - -24.29770992) / (7.587014215 - 7.229142741) 
    
  if( 7.587014215 < T_kK and T_kK < 8.241070825 ) : 
    log10Q_HI_neutral = -24.20610687 + (T_kK - 7.587014215) * (-24.17557252 - -24.20610687) / (8.241070825 - 7.587014215) 
    
  if( 8.241070825 < T_kK and T_kK < 8.954999469 ) : 
    log10Q_HI_neutral = -24.17557252 + (T_kK - 8.241070825) * (-24.11450382 - -24.17557252) / (8.954999469 - 8.241070825) 
    
  if( 8.954999469 < T_kK and T_kK < 9.610416807 ) : 
    log10Q_HI_neutral = -24.11450382 + (T_kK - 8.954999469) * (-23.99236641 - -24.11450382) / (9.610416807 - 8.954999469) 
    
  if( 9.610416807 < T_kK and T_kK < 10.38648536 ) : 
    log10Q_HI_neutral = -23.99236641 + (T_kK - 9.610416807) * (-23.7480916 - -23.99236641) / (10.38648536 - 9.610416807) 
    
  if( 10.38648536 < T_kK and T_kK < 11.34228342 ) : 
    log10Q_HI_neutral = -23.7480916 + (T_kK - 10.38648536) * (-23.40458015 - -23.7480916) / (11.34228342 - 10.38648536) 
    
  if( 11.34228342 < T_kK and T_kK < 11.94100375 ) : 
    log10Q_HI_neutral = -23.40458015 + (T_kK - 11.34228342) * (-23.09923664 - -23.40458015) / (11.94100375 - 11.34228342) 
    
  if( 11.94100375 < T_kK and T_kK < 12.89634823 ) : 
    log10Q_HI_neutral = -23.09923664 + (T_kK - 11.94100375) * (-22.78625954 - -23.09923664) / (12.89634823 - 11.94100375) 
    
  if( 12.89634823 < T_kK and T_kK < 13.73217543 ) : 
    log10Q_HI_neutral = -22.78625954 + (T_kK - 12.89634823) * (-22.51908397 - -22.78625954) / (13.73217543 - 12.89634823) 
    
  if( 13.73217543 < T_kK and T_kK < 14.8047693 ) : 
    log10Q_HI_neutral = -22.51908397 + (T_kK - 13.73217543) * (-22.3129771 - -22.51908397) / (14.8047693 - 13.73217543) 
    
  if( 14.8047693 < T_kK and T_kK < 15.93666824 ) : 
    log10Q_HI_neutral = -22.3129771 + (T_kK - 14.8047693) * (-22.11450382 - -22.3129771) / (15.93666824 - 14.8047693) 
    
  if( 15.93666824 < T_kK and T_kK < 16.9491633 ) : 
    log10Q_HI_neutral = -22.11450382 + (T_kK - 15.93666824) * (-21.95419847 - -22.11450382) / (16.9491633 - 15.93666824) 
    
  if( 16.9491633 < T_kK and T_kK < 17.8415741 ) : 
    log10Q_HI_neutral = -21.95419847 + (T_kK - 16.9491633) * (-21.8778626 - -21.95419847) / (17.8415741 - 16.9491633) 
    
  if( 17.8415741 < T_kK and T_kK < 18.9119001 ) : 
    log10Q_HI_neutral = -21.8778626 + (T_kK - 17.8415741) * (-21.82442748 - -21.8778626) / (18.9119001 - 17.8415741) 
    
  if( 18.9119001 < T_kK and T_kK < 20.21887938 ) : 
    log10Q_HI_neutral = -21.82442748 + (T_kK - 18.9119001) * (-21.83969466 - -21.82442748) / (20.21887938 - 18.9119001) 
    
  if( 20.21887938 < T_kK and T_kK < 22.11947626 ) : 
    log10Q_HI_neutral = -21.83969466 + (T_kK - 20.21887938) * (-21.89312977 - -21.83969466) / (22.11947626 - 20.21887938) 
    
  if( 22.11947626 < T_kK and T_kK < 23.30693827 ) : 
    log10Q_HI_neutral = -21.89312977 + (T_kK - 22.11947626) * (-21.95419847 - -21.89312977) / (23.30693827 - 22.11947626) 
    
  if( 23.30693827 < T_kK and T_kK < 24.61244342 ) : 
    log10Q_HI_neutral = -21.95419847 + (T_kK - 23.30693827) * (-22.06870229 - -21.95419847) / (24.61244342 - 23.30693827) 
    
  if( 24.61244342 < T_kK and T_kK < 25.79899827 ) : 
    log10Q_HI_neutral = -22.06870229 + (T_kK - 24.61244342) * (-22.19083969 - -22.06870229) / (25.79899827 - 24.61244342) 
    
  if( 25.79899827 < T_kK and T_kK < 27.46044721 ) : 
    log10Q_HI_neutral = -22.19083969 + (T_kK - 25.79899827) * (-22.34351145 - -22.19083969) / (27.46044721 - 25.79899827) 
    
  if( 27.46044721 < T_kK and T_kK < 29.7164209 ) : 
    log10Q_HI_neutral = -22.34351145 + (T_kK - 27.46044721) * (-22.47328244 - -22.34351145) / (29.7164209 - 27.46044721) 
    
  if( 29.7164209 < T_kK and T_kK < 32.50772768 ) : 
    log10Q_HI_neutral = -22.47328244 + (T_kK - 29.7164209) * (-22.5648855 - -22.47328244) / (32.50772768 - 29.7164209) 
    
  if( 32.50772768 < T_kK and T_kK < 34.88321865 ) : 
    log10Q_HI_neutral = -22.5648855 + (T_kK - 32.50772768) * (-22.64885496 - -22.5648855) / (34.88321865 - 32.50772768) 
    
  if( 34.88321865 < T_kK and T_kK < 37.31812808 ) : 
    log10Q_HI_neutral = -22.64885496 + (T_kK - 34.88321865) * (-22.73282443 - -22.64885496) / (37.31812808 - 34.88321865) 
    
  if( 37.31812808 < T_kK and T_kK < 38.80313595 ) : 
    log10Q_HI_neutral = -22.73282443 + (T_kK - 37.31812808) * (-22.76335878 - -22.73282443) / (38.80313595 - 37.31812808) 
    
  if( 38.80313595 < T_kK and T_kK < 41.47605939 ) : 
    log10Q_HI_neutral = -22.76335878 + (T_kK - 38.80313595) * (-22.82442748 - -22.76335878) / (41.47605939 - 38.80313595) 
    
  if( 41.47605939 < T_kK and T_kK < 44.20862808 ) : 
    log10Q_HI_neutral = -22.82442748 + (T_kK - 41.47605939) * (-22.87022901 - -22.82442748) / (44.20862808 - 41.47605939) 
    
  if( 44.20862808 < T_kK and T_kK < 46.70340954 ) : 
    log10Q_HI_neutral = -22.87022901 + (T_kK - 44.20862808) * (-22.92366412 - -22.87022901) / (46.70340954 - 44.20862808) 
    
  if( 46.70340954 < T_kK and T_kK < 49.97040417 ) : 
    log10Q_HI_neutral = -22.92366412 + (T_kK - 46.70340954) * (-22.99236641 - -22.92366412) / (49.97040417 - 46.70340954) 
    
  if( T_kK > 49.97040417 ) : 
    log10Q_HI_neutral = -22.99236641 
    
  if( T_kK < 2.070975919 ) : 
    log10Q_CaII = -21.59770115 
    
  if( 2.070975919 < T_kK and T_kK < 2.532319392 ) : 
    log10Q_CaII = -21.59770115 + (T_kK - 2.070975919) * (-20.98467433 - -21.59770115) / (2.532319392 - 2.070975919) 
    
  if( 2.532319392 < T_kK and T_kK < 2.887198986 ) : 
    log10Q_CaII = -20.98467433 + (T_kK - 2.532319392) * (-20.61685824 - -20.98467433) / (2.887198986 - 2.532319392) 
    
  if( 2.887198986 < T_kK and T_kK < 3.27756654 ) : 
    log10Q_CaII = -20.61685824 + (T_kK - 2.887198986) * (-20.24904215 - -20.61685824) / (3.27756654 - 2.887198986) 
    
  if( 3.27756654 < T_kK and T_kK < 3.951837769 ) : 
    log10Q_CaII = -20.24904215 + (T_kK - 3.27756654) * (-19.85057471 - -20.24904215) / (3.951837769 - 3.27756654) 
    
  if( 3.951837769 < T_kK and T_kK < 4.732572877 ) : 
    log10Q_CaII = -19.85057471 + (T_kK - 3.951837769) * (-19.51340996 - -19.85057471) / (4.732572877 - 3.951837769) 
    
  if( 4.732572877 < T_kK and T_kK < 5.761723701 ) : 
    log10Q_CaII = -19.51340996 + (T_kK - 4.732572877) * (-19.23754789 - -19.51340996) / (5.761723701 - 4.732572877) 
    
  if( 5.761723701 < T_kK and T_kK < 6.684410646 ) : 
    log10Q_CaII = -19.23754789 + (T_kK - 5.761723701) * (-19.03831418 - -19.23754789) / (6.684410646 - 5.761723701) 
    
  if( 6.684410646 < T_kK and T_kK < 7.961977186 ) : 
    log10Q_CaII = -19.03831418 + (T_kK - 6.684410646) * (-18.85440613 - -19.03831418) / (7.961977186 - 6.684410646) 
    
  if( 7.961977186 < T_kK and T_kK < 8.920152091 ) : 
    log10Q_CaII = -18.85440613 + (T_kK - 7.961977186) * (-18.73180077 - -18.85440613) / (8.920152091 - 7.961977186) 
    
  if( 8.920152091 < T_kK and T_kK < 10.16223067 ) : 
    log10Q_CaII = -18.73180077 + (T_kK - 8.920152091) * (-18.6091954 - -18.73180077) / (10.16223067 - 8.920152091) 
    
  if( 10.16223067 < T_kK and T_kK < 11.22686946 ) : 
    log10Q_CaII = -18.6091954 + (T_kK - 10.16223067) * (-18.51724138 - -18.6091954) / (11.22686946 - 10.16223067) 
    
  if( 11.22686946 < T_kK and T_kK < 12.50443599 ) : 
    log10Q_CaII = -18.51724138 + (T_kK - 11.22686946) * (-18.42528736 - -18.51724138) / (12.50443599 - 11.22686946) 
    
  if( 12.50443599 < T_kK and T_kK < 14.03041825 ) : 
    log10Q_CaII = -18.42528736 + (T_kK - 12.50443599) * (-18.348659 - -18.42528736) / (14.03041825 - 12.50443599) 
    
  if( 14.03041825 < T_kK and T_kK < 15.27249683 ) : 
    log10Q_CaII = -18.348659 + (T_kK - 14.03041825) * (-18.30268199 - -18.348659) / (15.27249683 - 14.03041825) 
    
  if( 15.27249683 < T_kK and T_kK < 16.97591888 ) : 
    log10Q_CaII = -18.30268199 + (T_kK - 15.27249683) * (-18.22605364 - -18.30268199) / (16.97591888 - 15.27249683) 
    
  if( 16.97591888 < T_kK and T_kK < 18.92775665 ) : 
    log10Q_CaII = -18.22605364 + (T_kK - 16.97591888) * (-18.16475096 - -18.22605364) / (18.92775665 - 16.97591888) 
    
  if( 18.92775665 < T_kK and T_kK < 21.66032953 ) : 
    log10Q_CaII = -18.16475096 + (T_kK - 18.92775665) * (-18.10344828 - -18.16475096) / (21.66032953 - 18.92775665) 
    
  if( 21.66032953 < T_kK and T_kK < 23.6121673 ) : 
    log10Q_CaII = -18.10344828 + (T_kK - 21.66032953) * (-18.07279693 - -18.10344828) / (23.6121673 - 21.66032953) 
    
  if( 23.6121673 < T_kK and T_kK < 26.66413181 ) : 
    log10Q_CaII = -18.07279693 + (T_kK - 23.6121673) * (-18.02681992 - -18.07279693) / (26.66413181 - 23.6121673) 
    
  if( 26.66413181 < T_kK and T_kK < 30.0 ) : 
    log10Q_CaII = -18.02681992 + (T_kK - 26.66413181) * (-17.99616858 - -18.02681992) / (30.0 - 26.66413181) 
    
  if( T_kK > 30.0 ) : 
    log10Q_CaII = -17.99616858 
    
  if( T_kK < 2.784781707 ) : 
    log10Q_MgII = -24.98467433 
    
  if( 2.784781707 < T_kK and T_kK < 2.965107255 ) : 
    log10Q_MgII = -24.98467433 + (T_kK - 2.784781707) * (-24.54022989 - -24.98467433) / (2.965107255 - 2.784781707) 
    
  if( 2.965107255 < T_kK and T_kK < 3.46898016 ) : 
    log10Q_MgII = -24.54022989 + (T_kK - 2.965107255) * (-23.52873563 - -24.54022989) / (3.46898016 - 2.965107255) 
    
  if( 3.46898016 < T_kK and T_kK < 3.79218586 ) : 
    log10Q_MgII = -23.52873563 + (T_kK - 3.46898016) * (-23.03831418 - -23.52873563) / (3.79218586 - 3.46898016) 
    
  if( 3.79218586 < T_kK and T_kK < 4.151197013 ) : 
    log10Q_MgII = -23.03831418 + (T_kK - 3.79218586) * (-22.51724138 - -23.03831418) / (4.151197013 - 3.79218586) 
    
  if( 4.151197013 < T_kK and T_kK < 5.118490861 ) : 
    log10Q_MgII = -22.51724138 + (T_kK - 4.151197013) * (-21.56704981 - -22.51724138) / (5.118490861 - 4.151197013) 
    
  if( 5.118490861 < T_kK and T_kK < 5.655367645 ) : 
    log10Q_MgII = -21.56704981 + (T_kK - 5.118490861) * (-21.1532567 - -21.56704981) / (5.655367645 - 5.118490861) 
    
  if( 5.655367645 < T_kK and T_kK < 6.334851258 ) : 
    log10Q_MgII = -21.1532567 + (T_kK - 5.655367645) * (-20.75478927 - -21.1532567) / (6.334851258 - 5.655367645) 
    
  if( 6.334851258 < T_kK and T_kK < 6.906781853 ) : 
    log10Q_MgII = -20.75478927 + (T_kK - 6.334851258) * (-20.4789272 - -20.75478927) / (6.906781853 - 6.334851258) 
    
  if( 6.906781853 < T_kK and T_kK < 7.871547453 ) : 
    log10Q_MgII = -20.4789272 + (T_kK - 6.906781853) * (-20.09578544 - -20.4789272) / (7.871547453 - 6.906781853) 
    
  if( 7.871547453 < T_kK and T_kK < 8.907308978 ) : 
    log10Q_MgII = -20.09578544 + (T_kK - 7.871547453) * (-19.78927203 - -20.09578544) / (8.907308978 - 7.871547453) 
    
  if( 8.907308978 < T_kK and T_kK < 9.907060058 ) : 
    log10Q_MgII = -19.78927203 + (T_kK - 8.907308978) * (-19.55938697 - -19.78927203) / (9.907060058 - 8.907308978) 
    
  if( 9.907060058 < T_kK and T_kK < 11.08515509 ) : 
    log10Q_MgII = -19.55938697 + (T_kK - 9.907060058) * (-19.32950192 - -19.55938697) / (11.08515509 - 9.907060058) 
    
  if( 11.08515509 < T_kK and T_kK < 12.4414574 ) : 
    log10Q_MgII = -19.32950192 + (T_kK - 11.08515509) * (-19.1302682 - -19.32950192) / (12.4414574 - 11.08515509) 
    
  if( 12.4414574 < T_kK and T_kK < 13.61914245 ) : 
    log10Q_MgII = -19.1302682 + (T_kK - 12.4414574) * (-18.99233716 - -19.1302682) / (13.61914245 - 12.4414574) 
    
  if( 13.61914245 < T_kK and T_kK < 15.58174586 ) : 
    log10Q_MgII = -18.99233716 + (T_kK - 13.61914245) * (-18.80842912 - -18.99233716) / (15.58174586 - 13.61914245) 
    
  if( 15.58174586 < T_kK and T_kK < 16.79482637 ) : 
    log10Q_MgII = -18.80842912 + (T_kK - 15.58174586) * (-18.73180077 - -18.80842912) / (16.79482637 - 15.58174586) 
    
  if( 16.79482637 < T_kK and T_kK < 18.3646631 ) : 
    log10Q_MgII = -18.73180077 + (T_kK - 16.79482637) * (-18.63984674 - -18.73180077) / (18.3646631 - 16.79482637) 
    
  if( 18.3646631 < T_kK and T_kK < 21.78961857 ) : 
    log10Q_MgII = -18.63984674 + (T_kK - 18.3646631) * (-18.47126437 - -18.63984674) / (21.78961857 - 18.3646631) 
    
  if( 21.78961857 < T_kK and T_kK < 25.14296313 ) : 
    log10Q_MgII = -18.47126437 + (T_kK - 21.78961857) * (-18.36398467 - -18.47126437) / (25.14296313 - 21.78961857) 
    
  if( 25.14296313 < T_kK and T_kK < 27.46177612 ) : 
    log10Q_MgII = -18.36398467 + (T_kK - 25.14296313) * (-18.28735632 - -18.36398467) / (27.46177612 - 25.14296313) 
    
  if( 27.46177612 < T_kK and T_kK < 29.17401469 ) : 
    log10Q_MgII = -18.28735632 + (T_kK - 27.46177612) * (-18.25670498 - -18.28735632) / (29.17401469 - 27.46177612) 
    
  if( 29.17401469 < T_kK and T_kK < 30.06580277 ) : 
    log10Q_MgII = -18.25670498 + (T_kK - 29.17401469) * (-18.24137931 - -18.25670498) / (30.06580277 - 29.17401469) 
    
  if( T_kK > 30.06580277 ) : 
    log10Q_MgII = -18.24137931 
    


  

  if( logtau < -1.894670995 ) : 
    esc_HI = 1.0 
    
  if( -1.894670995 < logtau and logtau < -1.746281156 ) : 
    esc_HI = 0.9933138489 + (logtau - -1.894670995) * (0.9932120667 - 0.9933138489) / (-1.746281156 - -1.894670995) 
    
  if( -1.746281156 < logtau and logtau < -1.597907861 ) : 
    esc_HI = 0.9932120667 + (logtau - -1.746281156) * (0.9915835513 - 0.9932120667) / (-1.597907861 - -1.746281156) 
    
  if( -1.597907861 < logtau and logtau < -1.449551107 ) : 
    esc_HI = 0.9915835513 + (logtau - -1.597907861) * (0.9884283028 - 0.9915835513) / (-1.449551107 - -1.597907861) 
    
  if( -1.449551107 < logtau and logtau < -1.301191046 ) : 
    esc_HI = 0.9884283028 + (logtau - -1.449551107) * (0.9855784009 - 0.9884283028) / (-1.301191046 - -1.449551107) 
    
  if( -1.301191046 < logtau and logtau < -1.152801207 ) : 
    esc_HI = 0.9855784009 + (logtau - -1.301191046) * (0.9854766187 - 0.9855784009) / (-1.152801207 - -1.301191046) 
    
  if( -1.152801207 < logtau and logtau < -1.004420191 ) : 
    esc_HI = 0.9854766187 + (logtau - -1.152801207) * (0.9845605788 - 0.9854766187) / (-1.004420191 - -1.152801207) 
    
  if( -1.004420191 < logtau and logtau < -0.8560524092 ) : 
    esc_HI = 0.9845605788 + (logtau - -1.004420191) * (0.9824231524 - 0.9845605788) / (-0.8560524092 - -1.004420191) 
    
  if( -0.8560524092 < logtau and logtau < -0.707695656 ) : 
    esc_HI = 0.9824231524 + (logtau - -0.8560524092) * (0.9792679039 - 0.9824231524) / (-0.707695656 - -0.8560524092) 
    
  if( -0.707695656 < logtau and logtau < -0.5593014059 ) : 
    esc_HI = 0.9792679039 + (logtau - -0.707695656) * (0.9795732505 - 0.9792679039) / (-0.5593014059 - -0.707695656) 
    
  if( -0.5593014059 < logtau and logtau < -0.410845396 ) : 
    esc_HI = 0.9795732505 + (logtau - -0.5593014059) * (0.9855784009 - 0.9795732505) / (-0.410845396 - -0.5593014059) 
    
  if( -0.410845396 < logtau and logtau < -0.2624015175 ) : 
    esc_HI = 0.9855784009 + (logtau - -0.410845396) * (0.990463947 - 0.9855784009) / (-0.2624015175 - -0.410845396) 
    
  if( -0.2624015175 < logtau and logtau < -0.1140524843 ) : 
    esc_HI = 0.990463947 + (logtau - -0.2624015175) * (0.986596223 - 0.990463947) / (-0.1140524843 - -0.2624015175) 
    
  if( -0.1140524843 < logtau and logtau < 0.03423920061 ) : 
    esc_HI = 0.986596223 + (logtau - -0.1140524843) * (0.9774358241 - 0.986596223) / (0.03423920061 - -0.1140524843) 
    
  if( 0.03423920061 < logtau and logtau < 0.1825187542 ) : 
    esc_HI = 0.9774358241 + (logtau - 0.03423920061) * (0.9671558208 - 0.9774358241) / (0.1825187542 - 0.03423920061) 
    
  if( 0.1825187542 < logtau and logtau < 0.3308049248 ) : 
    esc_HI = 0.9671558208 + (logtau - 0.1825187542) * (0.9574865108 - 0.9671558208) / (0.3308049248 - 0.1825187542) 
    
  if( 0.3308049248 < logtau and logtau < 0.4791307981 ) : 
    esc_HI = 0.9574865108 + (logtau - 0.3308049248) * (0.9514813604 - 0.9574865108) / (0.4791307981 - 0.3308049248) 
    
  if( 0.4791307981 < logtau and logtau < 0.6275018884 ) : 
    esc_HI = 0.9514813604 + (logtau - 0.4791307981) * (0.9496492806 - 0.9514813604) / (0.6275018884 - 0.4791307981) 
    
  if( 0.6275018884 < logtau and logtau < 0.7758873157 ) : 
    esc_HI = 0.9496492806 + (logtau - 0.6275018884) * (0.9491403695 - 0.9496492806) / (0.7758873157 - 0.6275018884) 
    
  if( 0.7758873157 < logtau and logtau < 0.9242330404 ) : 
    esc_HI = 0.9491403695 + (logtau - 0.7758873157) * (0.9449672989 - 0.9491403695) / (0.9242330404 - 0.7758873157) 
    
  if( 0.9242330404 < logtau and logtau < 1.072524725 ) : 
    esc_HI = 0.9449672989 + (logtau - 0.9242330404) * (0.9358068999 - 0.9449672989) / (1.072524725 - 0.9242330404) 
    
  if( 1.072524725 < logtau and logtau < 1.220781119 ) : 
    esc_HI = 0.9358068999 + (logtau - 1.072524725) * (0.9233894702 - 0.9358068999) / (1.220781119 - 1.072524725) 
    
  if( 1.220781119 < logtau and logtau < 1.369073907 ) : 
    esc_HI = 0.9233894702 + (logtau - 1.220781119) * (0.9143308535 - 0.9233894702) / (1.369073907 - 1.220781119) 
    
  if( 1.369073907 < logtau and logtau < 1.51742294 ) : 
    esc_HI = 0.9143308535 + (logtau - 1.369073907) * (0.9104631295 - 0.9143308535) / (1.51742294 - 1.369073907) 
    
  if( 1.51742294 < logtau and logtau < 1.665785207 ) : 
    esc_HI = 0.9104631295 + (logtau - 1.51742294) * (0.907816792 - 0.9104631295) / (1.665785207 - 1.51742294) 
    
  if( 1.665785207 < logtau and logtau < 1.814129829 ) : 
    esc_HI = 0.907816792 + (logtau - 1.665785207) * (0.9035419392 - 0.907816792) / (1.814129829 - 1.665785207) 
    
  if( 1.814129829 < logtau and logtau < 1.962476657 ) : 
    esc_HI = 0.9035419392 + (logtau - 1.814129829) * (0.8994706508 - 0.9035419392) / (1.962476657 - 1.814129829) 
    
  if( 1.962476657 < logtau and logtau < 2.110799222 ) : 
    esc_HI = 0.8994706508 + (logtau - 1.962476657) * (0.8931601537 - 0.8994706508) / (2.110799222 - 1.962476657) 
    
  if( 2.110799222 < logtau and logtau < 2.259028044 ) : 
    esc_HI = 0.8931601537 + (logtau - 2.110799222) * (0.8781981687 - 0.8931601537) / (2.259028044 - 2.110799222) 
    
  if( 2.259028044 < logtau and logtau < 2.393686763 ) : 
    esc_HI = 0.8781981687 + (logtau - 2.259028044) * (0.8558513189 - 0.8781981687) / (2.393686763 - 2.259028044) 
    
  if( 2.393686763 < logtau and logtau < 2.521579055 ) : 
    esc_HI = 0.8558513189 + (logtau - 2.393686763) * (0.8315310252 - 0.8558513189) / (2.521579055 - 2.393686763) 
    
  if( 2.521579055 < logtau and logtau < 2.662946484 ) : 
    esc_HI = 0.8315310252 + (logtau - 2.521579055) * (0.805831017 - 0.8315310252) / (2.662946484 - 2.521579055) 
    
  if( 2.662946484 < logtau and logtau < 2.804334646 ) : 
    esc_HI = 0.805831017 + (logtau - 2.662946484) * (0.7820445144 - 0.805831017) / (2.804334646 - 2.662946484) 
    
  if( 2.804334646 < logtau and logtau < 2.945740895 ) : 
    esc_HI = 0.7820445144 + (logtau - 2.804334646) * (0.75992724 - 0.7820445144) / (2.945740895 - 2.804334646) 
    
  if( 2.945740895 < logtau and logtau < 3.093938838 ) : 
    esc_HI = 0.75992724 + (logtau - 2.945740895) * (0.7421153532 - 0.75992724) / (3.093938838 - 2.945740895) 
    
  if( 3.093938838 < logtau and logtau < 3.242169866 ) : 
    esc_HI = 0.7421153532 + (logtau - 3.093938838) * (0.7273569326 - 0.7421153532) / (3.242169866 - 3.093938838) 
    
  if( 3.242169866 < logtau and logtau < 3.39038766 ) : 
    esc_HI = 0.7273569326 + (logtau - 3.242169866) * (0.7113771256 - 0.7273569326) / (3.39038766 - 3.242169866) 
    
  if( 3.39038766 < logtau and logtau < 3.52505177 ) : 
    esc_HI = 0.7113771256 + (logtau - 3.39038766) * (0.6895278777 - 0.7113771256) / (3.52505177 - 3.39038766) 
    
  if( 3.52505177 < logtau and logtau < 3.632697477 ) : 
    esc_HI = 0.6895278777 + (logtau - 3.52505177) * (0.6641501799 - 0.6895278777) / (3.632697477 - 3.52505177) 
    
  if( 3.632697477 < logtau and logtau < 3.720106034 ) : 
    esc_HI = 0.6641501799 + (logtau - 3.632697477) * (0.6385858813 - 0.6641501799) / (3.720106034 - 3.632697477) 
    
  if( 3.720106034 < logtau and logtau < 3.794035019 ) : 
    esc_HI = 0.6385858813 + (logtau - 3.720106034) * (0.6139919065 - 0.6385858813) / (3.794035019 - 3.720106034) 
    
  if( 3.794035019 < logtau and logtau < 3.854495351 ) : 
    esc_HI = 0.6139919065 + (logtau - 3.794035019) * (0.5913758993 - 0.6139919065) / (3.854495351 - 3.794035019) 
    
  if( 3.854495351 < logtau and logtau < 3.908176671 ) : 
    esc_HI = 0.5913758993 + (logtau - 3.854495351) * (0.565625 - 0.5913758993) / (3.908176671 - 3.854495351) 
    
  if( 3.908176671 < logtau and logtau < 3.961842828 ) : 
    esc_HI = 0.565625 + (logtau - 3.908176671) * (0.5384745953 - 0.565625) / (3.961842828 - 3.908176671) 
    
  if( 3.961842828 < logtau and logtau < 4.01551505 ) : 
    esc_HI = 0.5384745953 + (logtau - 3.961842828) * (0.5118839928 - 0.5384745953) / (4.01551505 - 3.961842828) 
    
  if( 4.01551505 < logtau and logtau < 4.069196371 ) : 
    esc_HI = 0.5118839928 + (logtau - 4.01551505) * (0.4861330935 - 0.5118839928) / (4.069196371 - 4.01551505) 
    
  if( 4.069196371 < logtau and logtau < 4.122901955 ) : 
    esc_HI = 0.4861330935 + (logtau - 4.069196371) * (0.4626214029 - 0.4861330935) / (4.122901955 - 4.069196371) 
    
  if( 4.122901955 < logtau and logtau < 4.190065274 ) : 
    esc_HI = 0.4626214029 + (logtau - 4.122901955) * (0.4361241007 - 0.4626214029) / (4.190065274 - 4.122901955) 
    
  if( 4.190065274 < logtau and logtau < 4.270732831 ) : 
    esc_HI = 0.4361241007 + (logtau - 4.190065274) * (0.4109330036 - 0.4361241007) / (4.270732831 - 4.190065274) 
    
  if( 4.270732831 < logtau and logtau < 4.358151786 ) : 
    esc_HI = 0.4109330036 + (logtau - 4.270732831) * (0.3863283659 - 0.4109330036) / (4.358151786 - 4.270732831) 
    
  if( 4.358151786 < logtau and logtau < 4.459048623 ) : 
    esc_HI = 0.3863283659 + (logtau - 4.358151786) * (0.3605974595 - 0.3863283659) / (4.459048623 - 4.358151786) 
    
  if( 4.459048623 < logtau and logtau < 4.566687254 ) : 
    esc_HI = 0.3605974595 + (logtau - 4.459048623) * (0.3345666592 - 0.3605974595) / (4.566687254 - 4.459048623) 
    
  if( 4.566687254 < logtau and logtau < 4.687825978 ) : 
    esc_HI = 0.3345666592 + (logtau - 4.566687254) * (0.3094595324 - 0.3345666592) / (4.687825978 - 4.566687254) 
    
  if( 4.687825978 < logtau and logtau < 4.82927976 ) : 
    esc_HI = 0.3094595324 + (logtau - 4.687825978) * (0.2917290713 - 0.3094595324) / (4.82927976 - 4.687825978) 
    
  if( 4.82927976 < logtau and logtau < 4.977508583 ) : 
    esc_HI = 0.2917290713 + (logtau - 4.82927976) * (0.2767670863 - 0.2917290713) / (4.977508583 - 4.82927976) 
    
  if( 4.977508583 < logtau and logtau < 5.125687777 ) : 
    esc_HI = 0.2767670863 + (logtau - 4.977508583) * (0.2572249019 - 0.2767670863) / (5.125687777 - 4.977508583) 
    
  if( 5.125687777 < logtau and logtau < 5.267075939 ) : 
    esc_HI = 0.2572249019 + (logtau - 5.125687777) * (0.2334383993 - 0.2572249019) / (5.267075939 - 5.125687777) 
    
  if( 5.267075939 < logtau and logtau < 5.394954483 ) : 
    esc_HI = 0.2334383993 + (logtau - 5.267075939) * (0.2078492206 - 0.2334383993) / (5.394954483 - 5.267075939) 
    
  if( 5.394954483 < logtau and logtau < 5.50936006 ) : 
    esc_HI = 0.2078492206 + (logtau - 5.394954483) * (0.1838399281 - 0.2078492206) / (5.50936006 - 5.394954483) 
    
  if( 5.50936006 < logtau and logtau < 5.610272278 ) : 
    esc_HI = 0.1838399281 + (logtau - 5.50936006) * (0.15952852 - 0.1838399281) / (5.610272278 - 5.50936006) 
    
  if( 5.610272278 < logtau and logtau < 5.704439453 ) : 
    esc_HI = 0.15952852 + (logtau - 5.610272278) * (0.135217112 - 0.15952852) / (5.704439453 - 5.610272278) 
    
  if( 5.704439453 < logtau and logtau < 5.805350371 ) : 
    esc_HI = 0.135217112 + (logtau - 5.704439453) * (0.1107857464 - 0.135217112) / (5.805350371 - 5.704439453) 
    
  if( 5.805350371 < logtau and logtau < 5.912996584 ) : 
    esc_HI = 0.1107857464 + (logtau - 5.805350371) * (0.08545469874 - 0.1107857464) / (5.912996584 - 5.805350371) 
    
  if( 5.912996584 < logtau and logtau < 6.013907285 ) : 
    esc_HI = 0.08545469874 + (logtau - 5.912996584) * (0.06100334018 - 0.08545469874) / (6.013907285 - 5.912996584) 
    
  if( 6.013907285 < logtau and logtau < 6.114821669 ) : 
    esc_HI = 0.06100334018 + (logtau - 6.013907285) * (0.03689186151 - 0.06100334018) / (6.114821669 - 6.013907285) 
    
  if( 6.114821669 < logtau and logtau < 6.242743402 ) : 
    esc_HI = 0.03689186151 + (logtau - 6.114821669) * (0.01528858731 - 0.03689186151) / (6.242743402 - 6.114821669) 
    
  if( 6.242743402 < logtau and logtau < 6.391035087 ) : 
    esc_HI = 0.01528858731 + (logtau - 6.242743402) * (0.006128188358 - 0.01528858731) / (6.391035087 - 6.242743402) 
    
  if( 6.391035087 < logtau and logtau < 6.539421617 ) : 
    esc_HI = 0.006128188358 + (logtau - 6.391035087) * (0.005721059516 - 0.006128188358) / (6.539421617 - 6.391035087) 
    
  if( 6.539421617 < logtau and logtau < 6.687813662 ) : 
    esc_HI = 0.005721059516 + (logtau - 6.539421617) * (0.005822841727 - 0.005721059516) / (6.687813662 - 6.539421617) 
    
  if( 6.687813662 < logtau and logtau < 6.836204603 ) : 
    esc_HI = 0.005822841727 + (logtau - 6.687813662) * (0.005822841727 - 0.005822841727) / (6.836204603 - 6.687813662) 
    
  if( 6.836204603 < logtau and logtau < 6.984593339 ) : 
    esc_HI = 0.005822841727 + (logtau - 6.836204603) * (0.005619277305 - 0.005822841727) / (6.984593339 - 6.836204603) 
    
  if( 6.984593339 < logtau and logtau < 7.132986487 ) : 
    esc_HI = 0.005619277305 + (logtau - 6.984593339) * (0.005822841727 - 0.005619277305) / (7.132986487 - 6.984593339) 
    
  if( 7.132986487 < logtau and logtau < 7.281377428 ) : 
    esc_HI = 0.005822841727 + (logtau - 7.132986487) * (0.005822841727 - 0.005822841727) / (7.281377428 - 7.132986487) 
    
  if( 7.281377428 < logtau and logtau < 7.429766164 ) : 
    esc_HI = 0.005822841727 + (logtau - 7.281377428) * (0.005619277305 - 0.005822841727) / (7.429766164 - 7.281377428) 
    
  if( 7.429766164 < logtau and logtau < 7.578159311 ) : 
    esc_HI = 0.005619277305 + (logtau - 7.429766164) * (0.005822841727 - 0.005619277305) / (7.578159311 - 7.429766164) 
    
  if( 7.578159311 < logtau and logtau < 7.726550253 ) : 
    esc_HI = 0.005822841727 + (logtau - 7.578159311) * (0.005822841727 - 0.005822841727) / (7.726550253 - 7.578159311) 
    
  if( 7.726550253 < logtau and logtau < 7.874941195 ) : 
    esc_HI = 0.005822841727 + (logtau - 7.726550253) * (0.005822841727 - 0.005822841727) / (7.874941195 - 7.726550253) 
    
  if( 7.874941195 < logtau and logtau < 7.962620685 ) : 
    esc_HI = 0.005822841727 + (logtau - 7.874941195) * (0.005263039568 - 0.005822841727) / (7.962620685 - 7.874941195) 
    
  if( logtau > 7.962620685 ) : 
    esc_HI = 0.0 
   


  logtau = math.log10(mc)
 
  if( logtau < -9.908828775 ) : 
    esc_CaII = 1.0 
    
  if( -9.908828775 < logtau and logtau < -9.745374031 ) : 
    esc_CaII = 0.9960242287 + (logtau - -9.908828775) * (0.9968666334 - 0.9960242287) / (-9.745374031 - -9.908828775) 
    
  if( -9.745374031 < logtau and logtau < -9.581919287 ) : 
    esc_CaII = 0.9968666334 + (logtau - -9.745374031) * (0.9972009422 - 0.9968666334) / (-9.581919287 - -9.745374031) 
    
  if( -9.581919287 < logtau and logtau < -9.418464543 ) : 
    esc_CaII = 0.9972009422 + (logtau - -9.581919287) * (0.9972303935 - 0.9972009422) / (-9.418464543 - -9.581919287) 
    
  if( -9.418464543 < logtau and logtau < -9.255009799 ) : 
    esc_CaII = 0.9972303935 + (logtau - -9.418464543) * (0.9972598449 - 0.9972303935) / (-9.255009799 - -9.418464543) 
    
  if( -9.255009799 < logtau and logtau < -9.091555055 ) : 
    esc_CaII = 0.9972598449 + (logtau - -9.255009799) * (0.9972892962 - 0.9972598449) / (-9.091555055 - -9.255009799) 
    
  if( -9.091555055 < logtau and logtau < -8.928100311 ) : 
    esc_CaII = 0.9972892962 + (logtau - -9.091555055) * (0.9973187475 - 0.9972892962) / (-8.928100311 - -9.091555055) 
    
  if( -8.928100311 < logtau and logtau < -8.764645567 ) : 
    esc_CaII = 0.9973187475 + (logtau - -8.928100311) * (0.9973481988 - 0.9973187475) / (-8.764645567 - -8.928100311) 
    
  if( -8.764645567 < logtau and logtau < -8.601190823 ) : 
    esc_CaII = 0.9973481988 + (logtau - -8.764645567) * (0.9973776501 - 0.9973481988) / (-8.601190823 - -8.764645567) 
    
  if( -8.601190823 < logtau and logtau < -8.437736079 ) : 
    esc_CaII = 0.9973776501 + (logtau - -8.601190823) * (0.9974071014 - 0.9973776501) / (-8.437736079 - -8.601190823) 
    
  if( -8.437736079 < logtau and logtau < -8.274281335 ) : 
    esc_CaII = 0.9974071014 + (logtau - -8.437736079) * (0.9974365527 - 0.9974071014) / (-8.274281335 - -8.437736079) 
    
  if( -8.274281335 < logtau and logtau < -8.110826591 ) : 
    esc_CaII = 0.9974365527 + (logtau - -8.274281335) * (0.997466004 - 0.9974365527) / (-8.110826591 - -8.274281335) 
    
  if( -8.110826591 < logtau and logtau < -7.947371847 ) : 
    esc_CaII = 0.997466004 + (logtau - -8.110826591) * (0.9974954553 - 0.997466004) / (-7.947371847 - -8.110826591) 
    
  if( -7.947371847 < logtau and logtau < -7.783917104 ) : 
    esc_CaII = 0.9974954553 + (logtau - -7.947371847) * (0.9975249066 - 0.9974954553) / (-7.783917104 - -7.947371847) 
    
  if( -7.783917104 < logtau and logtau < -7.62046236 ) : 
    esc_CaII = 0.9975249066 + (logtau - -7.783917104) * (0.9975543579 - 0.9975249066) / (-7.62046236 - -7.783917104) 
    
  if( -7.62046236 < logtau and logtau < -7.457007616 ) : 
    esc_CaII = 0.9975543579 + (logtau - -7.62046236) * (0.9975838092 - 0.9975543579) / (-7.457007616 - -7.62046236) 
    
  if( -7.457007616 < logtau and logtau < -7.293552872 ) : 
    esc_CaII = 0.9975838092 + (logtau - -7.457007616) * (0.9976132605 - 0.9975838092) / (-7.293552872 - -7.457007616) 
    
  if( -7.293552872 < logtau and logtau < -7.130098128 ) : 
    esc_CaII = 0.9976132605 + (logtau - -7.293552872) * (0.9972362351 - 0.9976132605) / (-7.130098128 - -7.293552872) 
    
  if( -7.130098128 < logtau and logtau < -6.966643384 ) : 
    esc_CaII = 0.9972362351 + (logtau - -7.130098128) * (0.9965543522 - 0.9972362351) / (-6.966643384 - -7.130098128) 
    
  if( -6.966643384 < logtau and logtau < -6.80318864 ) : 
    esc_CaII = 0.9965543522 + (logtau - -6.966643384) * (0.9965838035 - 0.9965543522) / (-6.80318864 - -6.966643384) 
    
  if( -6.80318864 < logtau and logtau < -6.639733896 ) : 
    esc_CaII = 0.9965838035 + (logtau - -6.80318864) * (0.9966132548 - 0.9965838035) / (-6.639733896 - -6.80318864) 
    
  if( -6.639733896 < logtau and logtau < -6.476279152 ) : 
    esc_CaII = 0.9966132548 + (logtau - -6.639733896) * (0.9891228873 - 0.9966132548) / (-6.476279152 - -6.639733896) 
    
  if( -6.476279152 < logtau and logtau < -6.312824408 ) : 
    esc_CaII = 0.9891228873 + (logtau - -6.476279152) * (0.9617151619 - 0.9891228873) / (-6.312824408 - -6.476279152) 
    
  if( -6.312824408 < logtau and logtau < -6.179088708 ) : 
    esc_CaII = 0.9617151619 + (logtau - -6.312824408) * (0.9417638325 - 0.9617151619) / (-6.179088708 - -6.312824408) 
    
  if( -6.179088708 < logtau and logtau < -6.089931575 ) : 
    esc_CaII = 0.9417638325 + (logtau - -6.179088708) * (0.9131000628 - 0.9417638325) / (-6.089931575 - -6.179088708) 
    
  if( -6.089931575 < logtau and logtau < -5.993344681 ) : 
    esc_CaII = 0.9131000628 + (logtau - -6.089931575) * (0.8828527357 - 0.9131000628) / (-5.993344681 - -6.089931575) 
    
  if( -5.993344681 < logtau and logtau < -5.85217922 ) : 
    esc_CaII = 0.8828527357 + (logtau - -5.993344681) * (0.8675971877 - 0.8828527357) / (-5.85217922 - -5.993344681) 
    
  if( -5.85217922 < logtau and logtau < -5.71101376 ) : 
    esc_CaII = 0.8675971877 + (logtau - -5.85217922) * (0.8485309208 - 0.8675971877) / (-5.71101376 - -5.85217922) 
    
  if( -5.71101376 < logtau and logtau < -5.599567343 ) : 
    esc_CaII = 0.8485309208 + (logtau - -5.71101376) * (0.8235399823 - 0.8485309208) / (-5.599567343 - -5.71101376) 
    
  if( -5.599567343 < logtau and logtau < -5.51041021 ) : 
    esc_CaII = 0.8235399823 + (logtau - -5.599567343) * (0.7998584555 - 0.8235399823) / (-5.51041021 - -5.599567343) 
    
  if( -5.51041021 < logtau and logtau < -5.436112599 ) : 
    esc_CaII = 0.7998584555 + (logtau - -5.51041021) * (0.7746093161 - 0.7998584555) / (-5.436112599 - -5.51041021) 
    
  if( -5.436112599 < logtau and logtau < -5.369244749 ) : 
    esc_CaII = 0.7746093161 + (logtau - -5.436112599) * (0.7492470569 - 0.7746093161) / (-5.369244749 - -5.436112599) 
    
  if( -5.369244749 < logtau and logtau < -5.3023769 ) : 
    esc_CaII = 0.7492470569 + (logtau - -5.369244749) * (0.7227669867 - 0.7492470569) / (-5.3023769 - -5.369244749) 
    
  if( -5.3023769 < logtau and logtau < -5.23550905 ) : 
    esc_CaII = 0.7227669867 + (logtau - -5.3023769) * (0.6953367773 - 0.7227669867) / (-5.23550905 - -5.3023769) 
    
  if( -5.23550905 < logtau and logtau < -5.176070961 ) : 
    esc_CaII = 0.6953367773 + (logtau - -5.23550905) * (0.6685200252 - 0.6953367773) / (-5.176070961 - -5.23550905) 
    
  if( -5.176070961 < logtau and logtau < -5.124062634 ) : 
    esc_CaII = 0.6685200252 + (logtau - -5.176070961) * (0.6425402926 - 0.6685200252) / (-5.124062634 - -5.176070961) 
    
  if( -5.124062634 < logtau and logtau < -5.072054306 ) : 
    esc_CaII = 0.6425402926 + (logtau - -5.124062634) * (0.6146043909 - 0.6425402926) / (-5.072054306 - -5.124062634) 
    
  if( -5.072054306 < logtau and logtau < -5.012616217 ) : 
    esc_CaII = 0.6146043909 + (logtau - -5.072054306) * (0.5841547534 - 0.6146043909) / (-5.012616217 - -5.072054306) 
    
  if( -5.012616217 < logtau and logtau < -4.953178128 ) : 
    esc_CaII = 0.5841547534 + (logtau - -5.012616217) * (0.5564996431 - 0.5841547534) / (-4.953178128 - -5.012616217) 
    
  if( -4.953178128 < logtau and logtau < -4.89374004 ) : 
    esc_CaII = 0.5564996431 + (logtau - -4.953178128) * (0.529682891 - 0.5564996431) / (-4.89374004 - -4.953178128) 
    
  if( -4.89374004 < logtau and logtau < -4.834301951 ) : 
    esc_CaII = 0.529682891 + (logtau - -4.89374004) * (0.5039839498 - 0.529682891) / (-4.834301951 - -4.89374004) 
    
  if( -4.834301951 < logtau and logtau < -4.767434101 ) : 
    esc_CaII = 0.5039839498 + (logtau - -4.834301951) * (0.4767214121 - 0.5039839498) / (-4.767434101 - -4.834301951) 
    
  if( -4.767434101 < logtau and logtau < -4.69313649 ) : 
    esc_CaII = 0.4767214121 + (logtau - -4.767434101) * (0.4516958348 - 0.4767214121) / (-4.69313649 - -4.767434101) 
    
  if( -4.69313649 < logtau and logtau < -4.618838879 ) : 
    esc_CaII = 0.4516958348 + (logtau - -4.69313649) * (0.4275645063 - 0.4516958348) / (-4.618838879 - -4.69313649) 
    
  if( -4.618838879 < logtau and logtau < -4.537111507 ) : 
    esc_CaII = 0.4275645063 + (logtau - -4.618838879) * (0.4017578001 - 0.4275645063) / (-4.537111507 - -4.618838879) 
    
  if( -4.537111507 < logtau and logtau < -4.447954374 ) : 
    esc_CaII = 0.4017578001 + (logtau - -4.537111507) * (0.3773683264 - 0.4017578001) / (-4.447954374 - -4.537111507) 
    
  if( -4.447954374 < logtau and logtau < -4.366227002 ) : 
    esc_CaII = 0.3773683264 + (logtau - -4.447954374) * (0.3502575075 - 0.3773683264) / (-4.366227002 - -4.447954374) 
    
  if( -4.366227002 < logtau and logtau < -4.291929392 ) : 
    esc_CaII = 0.3502575075 + (logtau - -4.366227002) * (0.3234434327 - 0.3502575075) / (-4.291929392 - -4.366227002) 
    
  if( -4.291929392 < logtau and logtau < -4.217631781 ) : 
    esc_CaII = 0.3234434327 + (logtau - -4.291929392) * (0.2957351093 - 0.3234434327) / (-4.217631781 - -4.291929392) 
    
  if( -4.217631781 < logtau and logtau < -4.14333417 ) : 
    esc_CaII = 0.2957351093 + (logtau - -4.217631781) * (0.2675796615 - 0.2957351093) / (-4.14333417 - -4.217631781) 
    
  if( -4.14333417 < logtau and logtau < -4.046747276 ) : 
    esc_CaII = 0.2675796615 + (logtau - -4.14333417) * (0.2401827522 - 0.2675796615) / (-4.046747276 - -4.14333417) 
    
  if( -4.046747276 < logtau and logtau < -3.905581815 ) : 
    esc_CaII = 0.2401827522 + (logtau - -4.046747276) * (0.2237839886 - 0.2401827522) / (-3.905581815 - -4.046747276) 
    
  if( -3.905581815 < logtau and logtau < -3.742127071 ) : 
    esc_CaII = 0.2237839886 + (logtau - -3.905581815) * (0.2111110433 - 0.2237839886) / (-3.742127071 - -3.905581815) 
    
  if( -3.742127071 < logtau and logtau < -3.578672327 ) : 
    esc_CaII = 0.2111110433 + (logtau - -3.742127071) * (0.1995559088 - 0.2111110433) / (-3.578672327 - -3.742127071) 
    
  if( -3.578672327 < logtau and logtau < -3.422647344 ) : 
    esc_CaII = 0.1995559088 + (logtau - -3.578672327) * (0.176049021 - 0.1995559088) / (-3.422647344 - -3.578672327) 
    
  if( -3.422647344 < logtau and logtau < -3.266622361 ) : 
    esc_CaII = 0.176049021 + (logtau - -3.422647344) * (0.1627243743 - 0.176049021) / (-3.266622361 - -3.422647344) 
    
  if( -3.266622361 < logtau and logtau < -3.103167617 ) : 
    esc_CaII = 0.1627243743 + (logtau - -3.266622361) * (0.164176494 - 0.1627243743) / (-3.103167617 - -3.266622361) 
    
  if( -3.103167617 < logtau and logtau < -2.939712873 ) : 
    esc_CaII = 0.164176494 + (logtau - -3.103167617) * (0.1641043261 - 0.164176494) / (-2.939712873 - -3.103167617) 
    
  if( -2.939712873 < logtau and logtau < -2.776258129 ) : 
    esc_CaII = 0.1641043261 + (logtau - -2.939712873) * (0.1590528188 - 0.1641043261) / (-2.776258129 - -2.939712873) 
    
  if( -2.776258129 < logtau and logtau < -2.649952191 ) : 
    esc_CaII = 0.1590528188 + (logtau - -2.776258129) * (0.1422067939 - 0.1590528188) / (-2.649952191 - -2.776258129) 
    
  if( -2.649952191 < logtau and logtau < -2.568224819 ) : 
    esc_CaII = 0.1422067939 + (logtau - -2.649952191) * (0.1153940579 - 0.1422067939) / (-2.568224819 - -2.649952191) 
    
  if( -2.568224819 < logtau and logtau < -2.449348641 ) : 
    esc_CaII = 0.1153940579 + (logtau - -2.568224819) * (0.09133173303 - 0.1153940579) / (-2.449348641 - -2.568224819) 
    
  if( -2.449348641 < logtau and logtau < -2.285893897 ) : 
    esc_CaII = 0.09133173303 + (logtau - -2.449348641) * (0.08008145614 - 0.09133173303) / (-2.285893897 - -2.449348641) 
    
  if( -2.285893897 < logtau and logtau < -2.122439153 ) : 
    esc_CaII = 0.08008145614 + (logtau - -2.285893897) * (0.05541744841 - 0.08008145614) / (-2.122439153 - -2.285893897) 
    
  if( -2.122439153 < logtau and logtau < -1.958984409 ) : 
    esc_CaII = 0.05541744841 + (logtau - -2.122439153) * (0.03573278016 - 0.05541744841) / (-1.958984409 - -2.122439153) 
    
  if( -1.958984409 < logtau and logtau < -1.795529665 ) : 
    esc_CaII = 0.03573278016 + (logtau - -1.958984409) * (0.02275497733 - 0.03573278016) / (-1.795529665 - -1.958984409) 
    
  if( -1.795529665 < logtau and logtau < -1.632074922 ) : 
    esc_CaII = 0.02275497733 + (logtau - -1.795529665) * (0.01556946735 - 0.02275497733) / (-1.632074922 - -1.795529665) 
    
  if( -1.632074922 < logtau and logtau < -1.468620178 ) : 
    esc_CaII = 0.01556946735 + (logtau - -1.632074922) * (0.01295682016 - 0.01556946735) / (-1.468620178 - -1.632074922) 
    
  if( -1.468620178 < logtau and logtau < -1.305165434 ) : 
    esc_CaII = 0.01295682016 + (logtau - -1.468620178) * (0.01166522222 - 0.01295682016) / (-1.305165434 - -1.468620178) 
    
  if( -1.305165434 < logtau and logtau < -1.14171069 ) : 
    esc_CaII = 0.01166522222 + (logtau - -1.305165434) * (0.01149143518 - 0.01166522222) / (-1.14171069 - -1.305165434) 
    
  if( -1.14171069 < logtau and logtau < -0.9782559456 ) : 
    esc_CaII = 0.01149143518 + (logtau - -1.14171069) * (0.01091117144 - 0.01149143518) / (-0.9782559456 - -1.14171069) 
    
  if( -0.9782559456 < logtau and logtau < -0.8148012017 ) : 
    esc_CaII = 0.01091117144 + (logtau - -0.9782559456) * (0.01053414606 - 0.01091117144) / (-0.8148012017 - -0.9782559456) 
    
  if( -0.8148012017 < logtau and logtau < -0.6513464577 ) : 
    esc_CaII = 0.01053414606 + (logtau - -0.8148012017) * (0.009649024806 - 0.01053414606) / (-0.6513464577 - -0.8148012017) 
    
  if( -0.6513464577 < logtau and logtau < -0.4878917137 ) : 
    esc_CaII = 0.009649024806 + (logtau - -0.6513464577) * (0.008459046035 - 0.009649024806) / (-0.4878917137 - -0.6513464577) 
    
  if( -0.4878917137 < logtau and logtau < -0.3244369698 ) : 
    esc_CaII = 0.008459046035 + (logtau - -0.4878917137) * (0.00848849734 - 0.008459046035) / (-0.3244369698 - -0.4878917137) 
    
  if( -0.3244369698 < logtau and logtau < -0.1609822258 ) : 
    esc_CaII = 0.00848849734 + (logtau - -0.3244369698) * (0.008517948646 - 0.00848849734) / (-0.1609822258 - -0.3244369698) 
    
  if( -0.1609822258 < logtau and logtau < 0.002472518157 ) : 
    esc_CaII = 0.008517948646 + (logtau - -0.1609822258) * (0.008445780778 - 0.008517948646) / (0.002472518157 - -0.1609822258) 
    
  if( 0.002472518157 < logtau and logtau < 0.1659272621 ) : 
    esc_CaII = 0.008445780778 + (logtau - 0.002472518157) * (0.008576851256 - 0.008445780778) / (0.1659272621 - 0.002472518157) 
    
  if( 0.1659272621 < logtau and logtau < 0.3293820061 ) : 
    esc_CaII = 0.008576851256 + (logtau - 0.1659272621) * (0.008606302561 - 0.008576851256) / (0.3293820061 - 0.1659272621) 
    
  if( 0.3293820061 < logtau and logtau < 0.4928367501 ) : 
    esc_CaII = 0.008606302561 + (logtau - 0.3293820061) * (0.008432515521 - 0.008606302561) / (0.4928367501 - 0.3293820061) 
    
  if( 0.4928367501 < logtau and logtau < 0.656291494 ) : 
    esc_CaII = 0.008432515521 + (logtau - 0.4928367501) * (0.008665205172 - 0.008432515521) / (0.656291494 - 0.4928367501) 
    
  if( 0.656291494 < logtau and logtau < 0.819746238 ) : 
    esc_CaII = 0.008665205172 + (logtau - 0.656291494) * (0.008084941439 - 0.008665205172) / (0.819746238 - 0.656291494) 
    
  if( 0.819746238 < logtau and logtau < 0.9386224154 ) : 
    esc_CaII = 0.008084941439 + (logtau - 0.819746238) * (0.008716075608 - 0.008084941439) / (0.9386224154 - 0.819746238) 
    
  if( logtau > 0.9386224154 ) : 
    esc_CaII = 0.0 
    
  if( logtau < -9.892532212 ) : 
    esc_MgII = 1.0 
    
  if( -9.892532212 < logtau and logtau < -9.729319419 ) : 
    esc_MgII = 0.9975049382 + (logtau - -9.892532212) * (0.9975049382 - 0.9975049382) / (-9.729319419 - -9.892532212) 
    
  if( -9.729319419 < logtau and logtau < -9.566106625 ) : 
    esc_MgII = 0.9975049382 + (logtau - -9.729319419) * (0.9975049382 - 0.9975049382) / (-9.566106625 - -9.729319419) 
    
  if( -9.566106625 < logtau and logtau < -9.402893832 ) : 
    esc_MgII = 0.9975049382 + (logtau - -9.566106625) * (0.9975049382 - 0.9975049382) / (-9.402893832 - -9.566106625) 
    
  if( -9.402893832 < logtau and logtau < -9.239681039 ) : 
    esc_MgII = 0.9975049382 + (logtau - -9.402893832) * (0.9975049382 - 0.9975049382) / (-9.239681039 - -9.402893832) 
    
  if( -9.239681039 < logtau and logtau < -9.076468246 ) : 
    esc_MgII = 0.9975049382 + (logtau - -9.239681039) * (0.9975049382 - 0.9975049382) / (-9.076468246 - -9.239681039) 
    
  if( -9.076468246 < logtau and logtau < -8.913255453 ) : 
    esc_MgII = 0.9975049382 + (logtau - -9.076468246) * (0.9975049382 - 0.9975049382) / (-8.913255453 - -9.076468246) 
    
  if( -8.913255453 < logtau and logtau < -8.75004266 ) : 
    esc_MgII = 0.9975049382 + (logtau - -8.913255453) * (0.9975049382 - 0.9975049382) / (-8.75004266 - -8.913255453) 
    
  if( -8.75004266 < logtau and logtau < -8.586829867 ) : 
    esc_MgII = 0.9975049382 + (logtau - -8.75004266) * (0.9975049382 - 0.9975049382) / (-8.586829867 - -8.75004266) 
    
  if( -8.586829867 < logtau and logtau < -8.423617074 ) : 
    esc_MgII = 0.9975049382 + (logtau - -8.586829867) * (0.9975049382 - 0.9975049382) / (-8.423617074 - -8.586829867) 
    
  if( -8.423617074 < logtau and logtau < -8.260404281 ) : 
    esc_MgII = 0.9975049382 + (logtau - -8.423617074) * (0.9975049382 - 0.9975049382) / (-8.260404281 - -8.423617074) 
    
  if( -8.260404281 < logtau and logtau < -8.097191488 ) : 
    esc_MgII = 0.9975049382 + (logtau - -8.260404281) * (0.9975049382 - 0.9975049382) / (-8.097191488 - -8.260404281) 
    
  if( -8.097191488 < logtau and logtau < -7.933978694 ) : 
    esc_MgII = 0.9975049382 + (logtau - -8.097191488) * (0.9975049382 - 0.9975049382) / (-7.933978694 - -8.097191488) 
    
  if( -7.933978694 < logtau and logtau < -7.770765901 ) : 
    esc_MgII = 0.9975049382 + (logtau - -7.933978694) * (0.9975049382 - 0.9975049382) / (-7.770765901 - -7.933978694) 
    
  if( -7.770765901 < logtau and logtau < -7.607553108 ) : 
    esc_MgII = 0.9975049382 + (logtau - -7.770765901) * (0.9975049382 - 0.9975049382) / (-7.607553108 - -7.770765901) 
    
  if( -7.607553108 < logtau and logtau < -7.444340315 ) : 
    esc_MgII = 0.9975049382 + (logtau - -7.607553108) * (0.9975049382 - 0.9975049382) / (-7.444340315 - -7.607553108) 
    
  if( -7.444340315 < logtau and logtau < -7.281127522 ) : 
    esc_MgII = 0.9975049382 + (logtau - -7.444340315) * (0.9975049382 - 0.9975049382) / (-7.281127522 - -7.444340315) 
    
  if( -7.281127522 < logtau and logtau < -7.117914729 ) : 
    esc_MgII = 0.9975049382 + (logtau - -7.281127522) * (0.9975049382 - 0.9975049382) / (-7.117914729 - -7.281127522) 
    
  if( -7.117914729 < logtau and logtau < -6.954701936 ) : 
    esc_MgII = 0.9975049382 + (logtau - -7.117914729) * (0.9975049382 - 0.9975049382) / (-6.954701936 - -7.117914729) 
    
  if( -6.954701936 < logtau and logtau < -6.791489143 ) : 
    esc_MgII = 0.9975049382 + (logtau - -6.954701936) * (0.9975049382 - 0.9975049382) / (-6.791489143 - -6.954701936) 
    
  if( -6.791489143 < logtau and logtau < -6.62827635 ) : 
    esc_MgII = 0.9975049382 + (logtau - -6.791489143) * (0.9975049382 - 0.9975049382) / (-6.62827635 - -6.791489143) 
    
  if( -6.62827635 < logtau and logtau < -6.465063557 ) : 
    esc_MgII = 0.9975049382 + (logtau - -6.62827635) * (0.9975049382 - 0.9975049382) / (-6.465063557 - -6.62827635) 
    
  if( -6.465063557 < logtau and logtau < -6.301850763 ) : 
    esc_MgII = 0.9975049382 + (logtau - -6.465063557) * (0.9975049382 - 0.9975049382) / (-6.301850763 - -6.465063557) 
    
  if( -6.301850763 < logtau and logtau < -6.13863797 ) : 
    esc_MgII = 0.9975049382 + (logtau - -6.301850763) * (0.9975049382 - 0.9975049382) / (-6.13863797 - -6.301850763) 
    
  if( -6.13863797 < logtau and logtau < -5.975425177 ) : 
    esc_MgII = 0.9975049382 + (logtau - -6.13863797) * (0.9975049382 - 0.9975049382) / (-5.975425177 - -6.13863797) 
    
  if( -5.975425177 < logtau and logtau < -5.812212384 ) : 
    esc_MgII = 0.9975049382 + (logtau - -5.975425177) * (0.9975049382 - 0.9975049382) / (-5.812212384 - -5.975425177) 
    
  if( -5.812212384 < logtau and logtau < -5.648999591 ) : 
    esc_MgII = 0.9975049382 + (logtau - -5.812212384) * (0.9975049382 - 0.9975049382) / (-5.648999591 - -5.812212384) 
    
  if( -5.648999591 < logtau and logtau < -5.485786798 ) : 
    esc_MgII = 0.9975049382 + (logtau - -5.648999591) * (0.9975049382 - 0.9975049382) / (-5.485786798 - -5.648999591) 
    
  if( -5.485786798 < logtau and logtau < -5.322574005 ) : 
    esc_MgII = 0.9975049382 + (logtau - -5.485786798) * (0.9975049382 - 0.9975049382) / (-5.322574005 - -5.485786798) 
    
  if( -5.322574005 < logtau and logtau < -5.159361212 ) : 
    esc_MgII = 0.9975049382 + (logtau - -5.322574005) * (0.9975049382 - 0.9975049382) / (-5.159361212 - -5.322574005) 
    
  if( -5.159361212 < logtau and logtau < -4.996148419 ) : 
    esc_MgII = 0.9975049382 + (logtau - -5.159361212) * (0.9975049382 - 0.9975049382) / (-4.996148419 - -5.159361212) 
    
  if( -4.996148419 < logtau and logtau < -4.832935626 ) : 
    esc_MgII = 0.9975049382 + (logtau - -4.996148419) * (0.9956719818 - 0.9975049382) / (-4.832935626 - -4.996148419) 
    
  if( -4.832935626 < logtau and logtau < -4.669722832 ) : 
    esc_MgII = 0.9956719818 + (logtau - -4.832935626) * (0.9855907215 - 0.9956719818) / (-4.669722832 - -4.832935626) 
    
  if( -4.669722832 < logtau and logtau < -4.528766329 ) : 
    esc_MgII = 0.9855907215 + (logtau - -4.669722832) * (0.9658609824 - 0.9855907215) / (-4.528766329 - -4.669722832) 
    
  if( -4.528766329 < logtau and logtau < -4.417484879 ) : 
    esc_MgII = 0.9658609824 + (logtau - -4.528766329) * (0.9424580568 - 0.9658609824) / (-4.417484879 - -4.528766329) 
    
  if( -4.417484879 < logtau and logtau < -4.321040956 ) : 
    esc_MgII = 0.9424580568 + (logtau - -4.417484879) * (0.9174149262 - 0.9424580568) / (-4.321040956 - -4.417484879) 
    
  if( -4.321040956 < logtau and logtau < -4.23943456 ) : 
    esc_MgII = 0.9174149262 + (logtau - -4.321040956) * (0.8926598316 - 0.9174149262) / (-4.23943456 - -4.321040956) 
    
  if( -4.23943456 < logtau and logtau < -4.17266569 ) : 
    esc_MgII = 0.8926598316 + (logtau - -4.23943456) * (0.870369045 - 0.8926598316) / (-4.17266569 - -4.23943456) 
    
  if( -4.17266569 < logtau and logtau < -4.113315583 ) : 
    esc_MgII = 0.870369045 + (logtau - -4.17266569) * (0.8448858594 - 0.870369045) / (-4.113315583 - -4.17266569) 
    
  if( -4.113315583 < logtau and logtau < -4.06138424 ) : 
    esc_MgII = 0.8448858594 + (logtau - -4.113315583) * (0.8223897139 - 0.8448858594) / (-4.06138424 - -4.113315583) 
    
  if( -4.06138424 < logtau and logtau < -4.009452897 ) : 
    esc_MgII = 0.8223897139 + (logtau - -4.06138424) * (0.7947595932 - 0.8223897139) / (-4.009452897 - -4.06138424) 
    
  if( -4.009452897 < logtau and logtau < -3.957521553 ) : 
    esc_MgII = 0.7947595932 + (logtau - -4.009452897) * (0.7671294725 - 0.7947595932) / (-3.957521553 - -4.009452897) 
    
  if( -3.957521553 < logtau and logtau < -3.913008973 ) : 
    esc_MgII = 0.7671294725 + (logtau - -3.957521553) * (0.7406194918 - 0.7671294725) / (-3.913008973 - -3.957521553) 
    
  if( -3.913008973 < logtau and logtau < -3.868496394 ) : 
    esc_MgII = 0.7406194918 + (logtau - -3.913008973) * (0.7118692311 - 0.7406194918) / (-3.868496394 - -3.913008973) 
    
  if( -3.868496394 < logtau and logtau < -3.823983814 ) : 
    esc_MgII = 0.7118692311 + (logtau - -3.868496394) * (0.6831189704 - 0.7118692311) / (-3.823983814 - -3.868496394) 
    
  if( -3.823983814 < logtau and logtau < -3.779471234 ) : 
    esc_MgII = 0.6831189704 + (logtau - -3.823983814) * (0.6588492698 - 0.6831189704) / (-3.779471234 - -3.823983814) 
    
  if( -3.779471234 < logtau and logtau < -3.72753989 ) : 
    esc_MgII = 0.6588492698 + (logtau - -3.779471234) * (0.6328993592 - 0.6588492698) / (-3.72753989 - -3.779471234) 
    
  if( -3.72753989 < logtau and logtau < -3.653352257 ) : 
    esc_MgII = 0.6328993592 + (logtau - -3.72753989) * (0.6052692385 - 0.6328993592) / (-3.653352257 - -3.72753989) 
    
  if( -3.653352257 < logtau and logtau < -3.527233281 ) : 
    esc_MgII = 0.6052692385 + (logtau - -3.653352257) * (0.5822384806 - 0.6052692385) / (-3.527233281 - -3.653352257) 
    
  if( -3.527233281 < logtau and logtau < -3.364020488 ) : 
    esc_MgII = 0.5822384806 + (logtau - -3.527233281) * (0.5803036933 - 0.5822384806) / (-3.364020488 - -3.527233281) 
    
  if( -3.364020488 < logtau and logtau < -3.200807694 ) : 
    esc_MgII = 0.5803036933 + (logtau - -3.364020488) * (0.5690004621 - 0.5803036933) / (-3.200807694 - -3.364020488) 
    
  if( -3.200807694 < logtau and logtau < -3.074688718 ) : 
    esc_MgII = 0.5690004621 + (logtau - -3.200807694) * (0.5432881569 - 0.5690004621) / (-3.074688718 - -3.200807694) 
    
  if( -3.074688718 < logtau and logtau < -2.993082321 ) : 
    esc_MgII = 0.5432881569 + (logtau - -3.074688718) * (0.5169648663 - 0.5432881569) / (-2.993082321 - -3.074688718) 
    
  if( -2.993082321 < logtau and logtau < -2.918894688 ) : 
    esc_MgII = 0.5169648663 + (logtau - -2.993082321) * (0.4916497017 - 0.5169648663) / (-2.918894688 - -2.993082321) 
    
  if( -2.918894688 < logtau and logtau < -2.837288292 ) : 
    esc_MgII = 0.4916497017 + (logtau - -2.918894688) * (0.464318285 - 0.4916497017) / (-2.837288292 - -2.918894688) 
    
  if( -2.837288292 < logtau and logtau < -2.748263132 ) : 
    esc_MgII = 0.464318285 + (logtau - -2.837288292) * (0.4378083043 - 0.464318285) / (-2.748263132 - -2.837288292) 
    
  if( -2.748263132 < logtau and logtau < -2.651819209 ) : 
    esc_MgII = 0.4378083043 + (logtau - -2.748263132) * (0.4161522638 - 0.4378083043) / (-2.651819209 - -2.748263132) 
    
  if( -2.651819209 < logtau and logtau < -2.577631575 ) : 
    esc_MgII = 0.4161522638 + (logtau - -2.651819209) * (0.384788343 - 0.4161522638) / (-2.577631575 - -2.651819209) 
    
  if( -2.577631575 < logtau and logtau < -2.540537759 ) : 
    esc_MgII = 0.384788343 + (logtau - -2.577631575) * (0.3562247723 - 0.384788343) / (-2.540537759 - -2.577631575) 
    
  if( -2.540537759 < logtau and logtau < -2.503443942 ) : 
    esc_MgII = 0.3562247723 + (logtau - -2.540537759) * (0.3272878216 - 0.3562247723) / (-2.503443942 - -2.540537759) 
    
  if( -2.503443942 < logtau and logtau < -2.458931362 ) : 
    esc_MgII = 0.3272878216 + (logtau - -2.503443942) * (0.2944303808 - 0.3272878216) / (-2.458931362 - -2.503443942) 
    
  if( -2.458931362 < logtau and logtau < -2.414418782 ) : 
    esc_MgII = 0.2944303808 + (logtau - -2.458931362) * (0.26343984 - 0.2944303808) / (-2.414418782 - -2.458931362) 
    
  if( -2.414418782 < logtau and logtau < -2.369906202 ) : 
    esc_MgII = 0.26343984 + (logtau - -2.414418782) * (0.2324492992 - 0.26343984) / (-2.369906202 - -2.414418782) 
    
  if( -2.369906202 < logtau and logtau < -2.325393622 ) : 
    esc_MgII = 0.2324492992 + (logtau - -2.369906202) * (0.2040724185 - 0.2324492992) / (-2.325393622 - -2.369906202) 
    
  if( -2.325393622 < logtau and logtau < -2.280881042 ) : 
    esc_MgII = 0.2040724185 + (logtau - -2.325393622) * (0.1753221578 - 0.2040724185) / (-2.280881042 - -2.325393622) 
    
  if( -2.280881042 < logtau and logtau < -2.177018356 ) : 
    esc_MgII = 0.1753221578 + (logtau - -2.280881042) * (0.1528175264 - 0.1753221578) / (-2.177018356 - -2.280881042) 
    
  if( -2.177018356 < logtau and logtau < -2.036061853 ) : 
    esc_MgII = 0.1528175264 + (logtau - -2.177018356) * (0.1377974669 - 0.1528175264) / (-2.036061853 - -2.177018356) 
    
  if( -2.036061853 < logtau and logtau < -1.932199166 ) : 
    esc_MgII = 0.1377974669 + (logtau - -2.036061853) * (0.1114741763 - 0.1377974669) / (-1.932199166 - -2.036061853) 
    
  if( -1.932199166 < logtau and logtau < -1.843174006 ) : 
    esc_MgII = 0.1114741763 + (logtau - -1.932199166) * (0.08552426561 - 0.1114741763) / (-1.843174006 - -1.932199166) 
    
  if( -1.843174006 < logtau and logtau < -1.731892557 ) : 
    esc_MgII = 0.08552426561 + (logtau - -1.843174006) * (0.06044557499 - 0.08552426561) / (-1.731892557 - -1.843174006) 
    
  if( -1.731892557 < logtau and logtau < -1.58351729 ) : 
    esc_MgII = 0.06044557499 + (logtau - -1.731892557) * (0.03774859625 - 0.06044557499) / (-1.58351729 - -1.731892557) 
    
  if( -1.58351729 < logtau and logtau < -1.420304497 ) : 
    esc_MgII = 0.03774859625 + (logtau - -1.58351729) * (0.01778973757 - 0.03774859625) / (-1.420304497 - -1.58351729) 
    
  if( -1.420304497 < logtau and logtau < -1.257091704 ) : 
    esc_MgII = 0.01778973757 + (logtau - -1.420304497) * (0.006079182735 - 0.01778973757) / (-1.257091704 - -1.420304497) 
    
  if( -1.257091704 < logtau and logtau < -1.093878911 ) : 
    esc_MgII = 0.006079182735 + (logtau - -1.257091704) * (0.005060873619 - 0.006079182735) / (-1.093878911 - -1.257091704) 
    
  if( -1.093878911 < logtau and logtau < -0.9306661177 ) : 
    esc_MgII = 0.005060873619 + (logtau - -1.093878911) * (0.005060873619 - 0.005060873619) / (-0.9306661177 - -1.093878911) 
    
  if( -0.9306661177 < logtau and logtau < -0.7674533246 ) : 
    esc_MgII = 0.005060873619 + (logtau - -0.9306661177) * (0.005060873619 - 0.005060873619) / (-0.7674533246 - -0.9306661177) 
    
  if( -0.7674533246 < logtau and logtau < -0.6042405315 ) : 
    esc_MgII = 0.005060873619 + (logtau - -0.7674533246) * (0.005060873619 - 0.005060873619) / (-0.6042405315 - -0.7674533246) 
    
  if( -0.6042405315 < logtau and logtau < -0.4410277384 ) : 
    esc_MgII = 0.005060873619 + (logtau - -0.6042405315) * (0.005060873619 - 0.005060873619) / (-0.4410277384 - -0.6042405315) 
    
  if( -0.4410277384 < logtau and logtau < -0.2778149453 ) : 
    esc_MgII = 0.005060873619 + (logtau - -0.4410277384) * (0.005060873619 - 0.005060873619) / (-0.2778149453 - -0.4410277384) 
    
  if( -0.2778149453 < logtau and logtau < -0.1146021522 ) : 
    esc_MgII = 0.005060873619 + (logtau - -0.2778149453) * (0.005060873619 - 0.005060873619) / (-0.1146021522 - -0.2778149453) 
    
  if( -0.1146021522 < logtau and logtau < 0.04861064087 ) : 
    esc_MgII = 0.005060873619 + (logtau - -0.1146021522) * (0.005060873619 - 0.005060873619) / (0.04861064087 - -0.1146021522) 
    
  if( 0.04861064087 < logtau and logtau < 0.211823434 ) : 
    esc_MgII = 0.005060873619 + (logtau - 0.04861064087) * (0.005060873619 - 0.005060873619) / (0.211823434 - 0.04861064087) 
    
  if( 0.211823434 < logtau and logtau < 0.3750362271 ) : 
    esc_MgII = 0.005060873619 + (logtau - 0.211823434) * (0.005060873619 - 0.005060873619) / (0.3750362271 - 0.211823434) 
    
  if( 0.3750362271 < logtau and logtau < 0.5382490202 ) : 
    esc_MgII = 0.005060873619 + (logtau - 0.3750362271) * (0.005060873619 - 0.005060873619) / (0.5382490202 - 0.3750362271) 
    
  if( 0.5382490202 < logtau and logtau < 0.7014618133 ) : 
    esc_MgII = 0.005060873619 + (logtau - 0.5382490202) * (0.005060873619 - 0.005060873619) / (0.7014618133 - 0.5382490202) 
    
  if( 0.7014618133 < logtau and logtau < 0.8646746064 ) : 
    esc_MgII = 0.005060873619 + (logtau - 0.7014618133) * (0.005060873619 - 0.005060873619) / (0.8646746064 - 0.7014618133) 
    
  if( 0.8646746064 < logtau and logtau < 0.9759560562 ) : 
    esc_MgII = 0.005060873619 + (logtau - 0.8646746064) * (0.005060873619 - 0.005060873619) / (0.9759560562 - 0.8646746064) 
    
  if( logtau > 0.9759560562 ) : 
    esc_MgII = 0.0 
    
  if( T_kK < 2.32844355 ) : 
    frac_HI = 1.0 
    
  if( 2.32844355 < T_kK and T_kK < 2.739035463 ) : 
    frac_HI = 1.04534562 + (T_kK - 2.32844355) * (1.045243865 - 1.04534562) / (2.739035463 - 2.32844355) 
    
  if( 2.739035463 < T_kK and T_kK < 3.149627377 ) : 
    frac_HI = 1.045243865 + (T_kK - 2.739035463) * (1.045243865 - 1.045243865) / (3.149627377 - 2.739035463) 
    
  if( 3.149627377 < T_kK and T_kK < 3.560219291 ) : 
    frac_HI = 1.045243865 + (T_kK - 3.149627377) * (1.045243865 - 1.045243865) / (3.560219291 - 3.149627377) 
    
  if( 3.560219291 < T_kK and T_kK < 3.970811204 ) : 
    frac_HI = 1.045243865 + (T_kK - 3.560219291) * (1.045243865 - 1.045243865) / (3.970811204 - 3.560219291) 
    
  if( 3.970811204 < T_kK and T_kK < 4.381403118 ) : 
    frac_HI = 1.045243865 + (T_kK - 3.970811204) * (1.045243865 - 1.045243865) / (4.381403118 - 3.970811204) 
    
  if( 4.381403118 < T_kK and T_kK < 4.791995032 ) : 
    frac_HI = 1.045243865 + (T_kK - 4.381403118) * (1.045243865 - 1.045243865) / (4.791995032 - 4.381403118) 
    
  if( 4.791995032 < T_kK and T_kK < 5.202586945 ) : 
    frac_HI = 1.045243865 + (T_kK - 4.791995032) * (1.045243865 - 1.045243865) / (5.202586945 - 4.791995032) 
    
  if( 5.202586945 < T_kK and T_kK < 5.613178859 ) : 
    frac_HI = 1.045243865 + (T_kK - 5.202586945) * (1.045243865 - 1.045243865) / (5.613178859 - 5.202586945) 
    
  if( 5.613178859 < T_kK and T_kK < 6.023770773 ) : 
    frac_HI = 1.045243865 + (T_kK - 5.613178859) * (1.041987699 - 1.045243865) / (6.023770773 - 5.613178859) 
    
  if( 6.023770773 < T_kK and T_kK < 6.359709611 ) : 
    frac_HI = 1.041987699 + (T_kK - 6.023770773) * (1.018860202 - 1.041987699) / (6.359709611 - 6.023770773) 
    
  if( 6.359709611 < T_kK and T_kK < 6.565005568 ) : 
    frac_HI = 1.018860202 + (T_kK - 6.359709611) * (0.9923566131 - 1.018860202) / (6.565005568 - 6.359709611) 
    
  if( 6.565005568 < T_kK and T_kK < 6.714311718 ) : 
    frac_HI = 0.9923566131 + (T_kK - 6.565005568) * (0.9638142866 - 0.9923566131) / (6.714311718 - 6.565005568) 
    
  if( 6.714311718 < T_kK and T_kK < 6.863617869 ) : 
    frac_HI = 0.9638142866 + (T_kK - 6.714311718) * (0.9361114403 - 0.9638142866) / (6.863617869 - 6.714311718) 
    
  if( 6.863617869 < T_kK and T_kK < 6.99426075 ) : 
    frac_HI = 0.9361114403 + (T_kK - 6.863617869) * (0.9094346254 - 0.9361114403) / (6.99426075 - 6.863617869) 
    
  if( 6.99426075 < T_kK and T_kK < 7.106240363 ) : 
    frac_HI = 0.9094346254 + (T_kK - 6.99426075) * (0.8855560777 - 0.9094346254) / (7.106240363 - 6.99426075) 
    
  if( 7.106240363 < T_kK and T_kK < 7.218219976 ) : 
    frac_HI = 0.8855560777 + (T_kK - 7.106240363) * (0.8616775301 - 0.8855560777) / (7.218219976 - 7.106240363) 
    
  if( 7.218219976 < T_kK and T_kK < 7.348862858 ) : 
    frac_HI = 0.8616775301 + (T_kK - 7.218219976) * (0.8362132976 - 0.8616775301) / (7.348862858 - 7.218219976) 
    
  if( 7.348862858 < T_kK and T_kK < 7.516832277 ) : 
    frac_HI = 0.8362132976 + (T_kK - 7.348862858) * (0.8133234711 - 0.8362132976) / (7.516832277 - 7.348862858) 
    
  if( 7.516832277 < T_kK and T_kK < 7.722128234 ) : 
    frac_HI = 0.8133234711 + (T_kK - 7.516832277) * (0.7864974153 - 0.8133234711) / (7.722128234 - 7.516832277) 
    
  if( 7.722128234 < T_kK and T_kK < 8.039403803 ) : 
    frac_HI = 0.7864974153 + (T_kK - 7.722128234) * (0.7620592142 - 0.7864974153) / (8.039403803 - 7.722128234) 
    
  if( 8.039403803 < T_kK and T_kK < 8.449995717 ) : 
    frac_HI = 0.7620592142 + (T_kK - 8.039403803) * (0.7415046689 - 0.7620592142) / (8.449995717 - 8.039403803) 
    
  if( 8.449995717 < T_kK and T_kK < 8.860587631 ) : 
    frac_HI = 0.7415046689 + (T_kK - 8.449995717) * (0.7263431479 - 0.7415046689) / (8.860587631 - 8.449995717) 
    
  if( 8.860587631 < T_kK and T_kK < 9.271179544 ) : 
    frac_HI = 0.7263431479 + (T_kK - 8.860587631) * (0.7094517889 - 0.7263431479) / (9.271179544 - 8.860587631) 
    
  if( 9.271179544 < T_kK and T_kK < 9.681771458 ) : 
    frac_HI = 0.7094517889 + (T_kK - 9.271179544) * (0.6929674506 - 0.7094517889) / (9.681771458 - 9.271179544) 
    
  if( 9.681771458 < T_kK and T_kK < 10.09236337 ) : 
    frac_HI = 0.6929674506 + (T_kK - 9.681771458) * (0.6763813571 - 0.6929674506) / (10.09236337 - 9.681771458) 
    
  if( 10.09236337 < T_kK and T_kK < 10.50295529 ) : 
    frac_HI = 0.6763813571 + (T_kK - 10.09236337) * (0.6592864878 - 0.6763813571) / (10.50295529 - 10.09236337) 
    
  if( 10.50295529 < T_kK and T_kK < 10.9135472 ) : 
    frac_HI = 0.6592864878 + (T_kK - 10.50295529) * (0.6429039047 - 0.6592864878) / (10.9135472 - 10.50295529) 
    
  if( 10.9135472 < T_kK and T_kK < 11.32413911 ) : 
    frac_HI = 0.6429039047 + (T_kK - 10.9135472) * (0.6267248319 - 0.6429039047) / (11.32413911 - 10.9135472) 
    
  if( 11.32413911 < T_kK and T_kK < 11.73473103 ) : 
    frac_HI = 0.6267248319 + (T_kK - 11.32413911) * (0.6094264523 - 0.6267248319) / (11.73473103 - 11.32413911) 
    
  if( 11.73473103 < T_kK and T_kK < 12.14532294 ) : 
    frac_HI = 0.6094264523 + (T_kK - 11.73473103) * (0.5917210519 - 0.6094264523) / (12.14532294 - 11.73473103) 
    
  if( 12.14532294 < T_kK and T_kK < 12.55591485 ) : 
    frac_HI = 0.5917210519 + (T_kK - 12.14532294) * (0.5732016101 - 0.5917210519) / (12.55591485 - 12.14532294) 
    
  if( 12.55591485 < T_kK and T_kK < 12.96650677 ) : 
    frac_HI = 0.5732016101 + (T_kK - 12.55591485) * (0.5526470648 - 0.5732016101) / (12.96650677 - 12.55591485) 
    
  if( 12.96650677 < T_kK and T_kK < 13.37709868 ) : 
    frac_HI = 0.5526470648 + (T_kK - 12.96650677) * (0.5286328436 - 0.5526470648) / (13.37709868 - 12.96650677) 
    
  if( 13.37709868 < T_kK and T_kK < 13.78769059 ) : 
    frac_HI = 0.5286328436 + (T_kK - 13.37709868) * (0.5050256431 - 0.5286328436) / (13.78769059 - 13.37709868) 
    
  if( 13.78769059 < T_kK and T_kK < 14.19828251 ) : 
    frac_HI = 0.5050256431 + (T_kK - 13.78769059) * (0.4811131771 - 0.5050256431) / (14.19828251 - 13.78769059) 
    
  if( 14.19828251 < T_kK and T_kK < 14.60887442 ) : 
    frac_HI = 0.4811131771 + (T_kK - 14.19828251) * (0.4575059766 - 0.4811131771) / (14.60887442 - 14.19828251) 
    
  if( 14.60887442 < T_kK and T_kK < 15.01946634 ) : 
    frac_HI = 0.4575059766 + (T_kK - 14.60887442) * (0.4340005313 - 0.4575059766) / (15.01946634 - 14.60887442) 
    
  if( 15.01946634 < T_kK and T_kK < 15.43005825 ) : 
    frac_HI = 0.4340005313 + (T_kK - 15.01946634) * (0.4103933307 - 0.4340005313) / (15.43005825 - 15.01946634) 
    
  if( 15.43005825 < T_kK and T_kK < 15.84065016 ) : 
    frac_HI = 0.4103933307 + (T_kK - 15.43005825) * (0.3903475613 - 0.4103933307) / (15.84065016 - 15.43005825) 
    
  if( 15.84065016 < T_kK and T_kK < 16.25124208 ) : 
    frac_HI = 0.3903475613 + (T_kK - 15.84065016) * (0.3705053023 - 0.3903475613) / (16.25124208 - 15.84065016) 
    
  if( 16.25124208 < T_kK and T_kK < 16.64317072 ) : 
    frac_HI = 0.3705053023 + (T_kK - 16.25124208) * (0.3467963466 - 0.3705053023) / (16.64317072 - 16.25124208) 
    
  if( 16.64317072 < T_kK and T_kK < 17.0164361 ) : 
    frac_HI = 0.3467963466 + (T_kK - 16.64317072) * (0.3218358023 - 0.3467963466) / (17.0164361 - 16.64317072) 
    
  if( 17.0164361 < T_kK and T_kK < 17.3710382 ) : 
    frac_HI = 0.3218358023 + (T_kK - 17.0164361) * (0.2980443119 - 0.3218358023) / (17.3710382 - 17.0164361) 
    
  if( 17.3710382 < T_kK and T_kK < 17.70697704 ) : 
    frac_HI = 0.2980443119 + (T_kK - 17.3710382) * (0.2730464573 - 0.2980443119) / (17.70697704 - 17.3710382) 
    
  if( 17.70697704 < T_kK and T_kK < 18.02425261 ) : 
    frac_HI = 0.2730464573 + (T_kK - 17.70697704) * (0.2488569911 - 0.2730464573) / (18.02425261 - 17.70697704) 
    
  if( 18.02425261 < T_kK and T_kK < 18.32286491 ) : 
    frac_HI = 0.2488569911 + (T_kK - 18.02425261) * (0.2259111992 - 0.2488569911) / (18.32286491 - 18.02425261) 
    
  if( 18.32286491 < T_kK and T_kK < 18.65880375 ) : 
    frac_HI = 0.2259111992 + (T_kK - 18.32286491) * (0.201622239 - 0.2259111992) / (18.65880375 - 18.32286491) 
    
  if( 18.65880375 < T_kK and T_kK < 19.03206913 ) : 
    frac_HI = 0.201622239 + (T_kK - 18.65880375) * (0.1783406551 - 0.201622239) / (19.03206913 - 18.65880375) 
    
  if( 19.03206913 < T_kK and T_kK < 19.42399777 ) : 
    frac_HI = 0.1783406551 + (T_kK - 19.03206913) * (0.1542755563 - 0.1783406551) / (19.42399777 - 19.03206913) 
    
  if( 19.42399777 < T_kK and T_kK < 19.83458969 ) : 
    frac_HI = 0.1542755563 + (T_kK - 19.42399777) * (0.131380642 - 0.1542755563) / (19.83458969 - 19.42399777) 
    
  if( 19.83458969 < T_kK and T_kK < 20.2451816 ) : 
    frac_HI = 0.131380642 + (T_kK - 19.83458969) * (0.1158121003 - 0.131380642) / (20.2451816 - 19.83458969) 
    
  if( 20.2451816 < T_kK and T_kK < 20.65577351 ) : 
    frac_HI = 0.1158121003 + (T_kK - 20.2451816) * (0.1006505793 - 0.1158121003) / (20.65577351 - 20.2451816) 
    
  if( 20.65577351 < T_kK and T_kK < 21.06636543 ) : 
    frac_HI = 0.1006505793 + (T_kK - 20.65577351) * (0.08660836517 - 0.1006505793) / (21.06636543 - 20.65577351) 
    
  if( 21.06636543 < T_kK and T_kK < 21.47695734 ) : 
    frac_HI = 0.08660836517 + (T_kK - 21.06636543) * (0.07643284771 - 0.08660836517) / (21.47695734 - 21.06636543) 
    
  if( 21.47695734 < T_kK and T_kK < 21.88754925 ) : 
    frac_HI = 0.07643284771 + (T_kK - 21.47695734) * (0.0668678613 - 0.07643284771) / (21.88754925 - 21.47695734) 
    
  if( 21.88754925 < T_kK and T_kK < 22.29814117 ) : 
    frac_HI = 0.0668678613 + (T_kK - 21.88754925) * (0.05730287489 - 0.0668678613) / (22.29814117 - 21.88754925) 
    
  if( 22.29814117 < T_kK and T_kK < 22.70873308 ) : 
    frac_HI = 0.05730287489 + (T_kK - 22.29814117) * (0.04987474714 - 0.05730287489) / (22.70873308 - 22.29814117) 
    
  if( 22.70873308 < T_kK and T_kK < 23.119325 ) : 
    frac_HI = 0.04987474714 + (T_kK - 22.70873308) * (0.04315890562 - 0.04987474714) / (23.119325 - 22.70873308) 
    
  if( 23.119325 < T_kK and T_kK < 23.52991691 ) : 
    frac_HI = 0.04315890562 + (T_kK - 23.119325) * (0.0368500848 - 0.04315890562) / (23.52991691 - 23.119325) 
    
  if( 23.52991691 < T_kK and T_kK < 23.94050882 ) : 
    frac_HI = 0.0368500848 + (T_kK - 23.52991691) * (0.03105003984 - 0.0368500848) / (23.94050882 - 23.52991691) 
    
  if( 23.94050882 < T_kK and T_kK < 24.35110074 ) : 
    frac_HI = 0.03105003984 + (T_kK - 23.94050882) * (0.02708158803 - 0.03105003984) / (24.35110074 - 23.94050882) 
    
  if( 24.35110074 < T_kK and T_kK < 24.76169265 ) : 
    frac_HI = 0.02708158803 + (T_kK - 24.35110074) * (0.02311313623 - 0.02708158803) / (24.76169265 - 24.35110074) 
    
  if( 24.76169265 < T_kK and T_kK < 25.17228456 ) : 
    frac_HI = 0.02311313623 + (T_kK - 24.76169265) * (0.01924643959 - 0.02311313623) / (25.17228456 - 24.76169265) 
    
  if( 25.17228456 < T_kK and T_kK < 25.58287648 ) : 
    frac_HI = 0.01924643959 + (T_kK - 25.17228456) * (0.01629553953 - 0.01924643959) / (25.58287648 - 25.17228456) 
    
  if( 25.58287648 < T_kK and T_kK < 25.99346839 ) : 
    frac_HI = 0.01629553953 + (T_kK - 25.58287648) * (0.01497272226 - 0.01629553953) / (25.99346839 - 25.58287648) 
    
  if( 25.99346839 < T_kK and T_kK < 26.4040603 ) : 
    frac_HI = 0.01497272226 + (T_kK - 25.99346839) * (0.01324288429 - 0.01497272226) / (26.4040603 - 25.99346839) 
    
  if( 26.4040603 < T_kK and T_kK < 26.81465222 ) : 
    frac_HI = 0.01324288429 + (T_kK - 26.4040603) * (0.01151304632 - 0.01324288429) / (26.81465222 - 26.4040603) 
    
  if( 26.81465222 < T_kK and T_kK < 27.22524413 ) : 
    frac_HI = 0.01151304632 + (T_kK - 26.81465222) * (0.01019022905 - 0.01151304632) / (27.22524413 - 26.81465222) 
    
  if( 27.22524413 < T_kK and T_kK < 27.63583605 ) : 
    frac_HI = 0.01019022905 + (T_kK - 27.22524413) * (0.009579698004 - 0.01019022905) / (27.63583605 - 27.22524413) 
    
  if( 27.63583605 < T_kK and T_kK < 28.04642796 ) : 
    frac_HI = 0.009579698004 + (T_kK - 27.63583605) * (0.008765656607 - 0.009579698004) / (28.04642796 - 27.63583605) 
    
  if( 28.04642796 < T_kK and T_kK < 28.45701987 ) : 
    frac_HI = 0.008765656607 + (T_kK - 28.04642796) * (0.007748104861 - 0.008765656607) / (28.45701987 - 28.04642796) 
    
  if( 28.45701987 < T_kK and T_kK < 28.86761179 ) : 
    frac_HI = 0.007748104861 + (T_kK - 28.45701987) * (0.006527042766 - 0.007748104861) / (28.86761179 - 28.45701987) 
    
  if( 28.86761179 < T_kK and T_kK < 29.2782037 ) : 
    frac_HI = 0.006527042766 + (T_kK - 28.86761179) * (0.006527042766 - 0.006527042766) / (29.2782037 - 28.86761179) 
    
  if( 29.2782037 < T_kK and T_kK < 29.68879561 ) : 
    frac_HI = 0.006527042766 + (T_kK - 29.2782037) * (0.005407735845 - 0.006527042766) / (29.68879561 - 29.2782037) 
    
  if( 29.68879561 < T_kK and T_kK < 29.95008138 ) : 
    frac_HI = 0.005407735845 + (T_kK - 29.68879561) * (0.005407735845 - 0.005407735845) / (29.95008138 - 29.68879561) 
    
  if( T_kK > 29.95008138 ) : 
    frac_HI = 0.0 
    
  if( T_kK < 2.350135199 ) : 
    frac_CaII = 1.0 
    
  if( 2.350135199 < T_kK and T_kK < 2.760694502 ) : 
    frac_CaII = 1.045138863 + (T_kK - 2.350135199) * (1.045240545 - 1.045138863) / (2.760694502 - 2.350135199) 
    
  if( 2.760694502 < T_kK and T_kK < 3.171250697 ) : 
    frac_CaII = 1.045240545 + (T_kK - 2.760694502) * (1.045240545 - 1.045240545) / (3.171250697 - 2.760694502) 
    
  if( 3.171250697 < T_kK and T_kK < 3.581806893 ) : 
    frac_CaII = 1.045240545 + (T_kK - 3.171250697) * (1.045240545 - 1.045240545) / (3.581806893 - 3.171250697) 
    
  if( 3.581806893 < T_kK and T_kK < 3.992363089 ) : 
    frac_CaII = 1.045240545 + (T_kK - 3.581806893) * (1.045240545 - 1.045240545) / (3.992363089 - 3.581806893) 
    
  if( 3.992363089 < T_kK and T_kK < 4.40291307 ) : 
    frac_CaII = 1.045240545 + (T_kK - 3.992363089) * (1.045037181 - 1.045240545) / (4.40291307 - 3.992363089) 
    
  if( 4.40291307 < T_kK and T_kK < 4.813400908 ) : 
    frac_CaII = 1.045037181 + (T_kK - 4.40291307) * (1.042800176 - 1.045037181) / (4.813400908 - 4.40291307) 
    
  if( 4.813400908 < T_kK and T_kK < 5.223885639 ) : 
    frac_CaII = 1.042800176 + (T_kK - 4.813400908) * (1.040461489 - 1.042800176) / (5.223885639 - 4.813400908) 
    
  if( 5.223885639 < T_kK and T_kK < 5.634376584 ) : 
    frac_CaII = 1.040461489 + (T_kK - 5.223885639) * (1.038326167 - 1.040461489) / (5.634376584 - 5.223885639) 
    
  if( 5.634376584 < T_kK and T_kK < 6.044861315 ) : 
    frac_CaII = 1.038326167 + (T_kK - 5.634376584) * (1.03598748 - 1.038326167) / (6.044861315 - 5.634376584) 
    
  if( 6.044861315 < T_kK and T_kK < 6.455342939 ) : 
    frac_CaII = 1.03598748 + (T_kK - 6.044861315) * (1.033547111 - 1.03598748) / (6.455342939 - 6.044861315) 
    
  if( 6.455342939 < T_kK and T_kK < 6.865824563 ) : 
    frac_CaII = 1.033547111 + (T_kK - 6.455342939) * (1.031106742 - 1.033547111) / (6.865824563 - 6.455342939) 
    
  if( 6.865824563 < T_kK and T_kK < 7.276315508 ) : 
    frac_CaII = 1.031106742 + (T_kK - 6.865824563) * (1.028971419 - 1.031106742) / (7.276315508 - 6.865824563) 
    
  if( 7.276315508 < T_kK and T_kK < 7.686800239 ) : 
    frac_CaII = 1.028971419 + (T_kK - 7.276315508) * (1.026632732 - 1.028971419) / (7.686800239 - 7.276315508) 
    
  if( 7.686800239 < T_kK and T_kK < 8.097309827 ) : 
    frac_CaII = 1.026632732 + (T_kK - 7.686800239) * (1.025107501 - 1.026632732) / (8.097309827 - 7.686800239) 
    
  if( 8.097309827 < T_kK and T_kK < 8.507788343 ) : 
    frac_CaII = 1.025107501 + (T_kK - 8.097309827) * (1.02256545 - 1.025107501) / (8.507788343 - 8.097309827) 
    
  if( 8.507788343 < T_kK and T_kK < 8.918186074 ) : 
    frac_CaII = 1.02256545 + (T_kK - 8.507788343) * (1.017379666 - 1.02256545) / (8.918186074 - 8.507788343) 
    
  if( 8.918186074 < T_kK and T_kK < 9.32845641 ) : 
    frac_CaII = 1.017379666 + (T_kK - 8.918186074) * (1.008024918 - 1.017379666) / (9.32845641 - 8.918186074) 
    
  if( 9.32845641 < T_kK and T_kK < 9.738518566 ) : 
    frac_CaII = 1.008024918 + (T_kK - 9.32845641) * (0.9918574736 - 1.008024918) / (9.738518566 - 9.32845641) 
    
  if( 9.738518566 < T_kK and T_kK < 10.12977832 ) : 
    frac_CaII = 0.9918574736 + (T_kK - 9.738518566) * (0.9710838326 - 0.9918574736) / (10.12977832 - 9.738518566) 
    
  if( 10.12977832 < T_kK and T_kK < 10.48358132 ) : 
    frac_CaII = 0.9710838326 + (T_kK - 10.12977832) * (0.9459423829 - 0.9710838326) / (10.48358132 - 10.12977832) 
    
  if( 10.48358132 < T_kK and T_kK < 10.76272368 ) : 
    frac_CaII = 0.9459423829 + (T_kK - 10.48358132) * (0.9203411045 - 0.9459423829) / (10.76272368 - 10.48358132) 
    
  if( 10.76272368 < T_kK and T_kK < 10.98586592 ) : 
    frac_CaII = 0.9203411045 + (T_kK - 10.76272368) * (0.8942427138 - 0.9203411045) / (10.98586592 - 10.76272368) 
    
  if( 10.98586592 < T_kK and T_kK < 11.20898537 ) : 
    frac_CaII = 0.8942427138 + (T_kK - 10.98586592) * (0.8673986549 - 0.8942427138) / (11.20898537 - 10.98586592) 
    
  if( 11.20898537 < T_kK and T_kK < 11.3948328 ) : 
    frac_CaII = 0.8673986549 + (T_kK - 11.20898537) * (0.8422323496 - 0.8673986549) / (11.3948328 - 11.20898537) 
    
  if( 11.3948328 < T_kK and T_kK < 11.56199124 ) : 
    frac_CaII = 0.8422323496 + (T_kK - 11.3948328) * (0.8161712424 - 0.8422323496) / (11.56199124 - 11.3948328) 
    
  if( 11.56199124 < T_kK and T_kK < 11.74772588 ) : 
    frac_CaII = 0.8161712424 + (T_kK - 11.56199124) * (0.787313879 - 0.8161712424) / (11.74772588 - 11.56199124) 
    
  if( 11.74772588 < T_kK and T_kK < 11.93345368 ) : 
    frac_CaII = 0.787313879 + (T_kK - 11.74772588) * (0.7582328151 - 0.787313879) / (11.93345368 - 11.74772588) 
    
  if( 11.93345368 < T_kK and T_kK < 12.10055743 ) : 
    frac_CaII = 0.7582328151 + (T_kK - 11.93345368) * (0.730382104 - 0.7582328151) / (12.10055743 - 11.93345368) 
    
  if( 12.10055743 < T_kK and T_kK < 12.24893631 ) : 
    frac_CaII = 0.730382104 + (T_kK - 12.10055743) * (0.7004621633 - 0.730382104) / (12.24893631 - 12.10055743) 
    
  if( 12.24893631 < T_kK and T_kK < 12.39734082 ) : 
    frac_CaII = 0.7004621633 + (T_kK - 12.24893631) * (0.6713810994 - 0.7004621633) / (12.39734082 - 12.24893631) 
    
  if( 12.39734082 < T_kK and T_kK < 12.54574534 ) : 
    frac_CaII = 0.6713810994 + (T_kK - 12.39734082) * (0.6423000355 - 0.6713810994) / (12.54574534 - 12.39734082) 
    
  if( 12.54574534 < T_kK and T_kK < 12.6941413 ) : 
    frac_CaII = 0.6423000355 + (T_kK - 12.54574534) * (0.612939346 - 0.6423000355) / (12.6941413 - 12.54574534) 
    
  if( 12.6941413 < T_kK and T_kK < 12.84253727 ) : 
    frac_CaII = 0.612939346 + (T_kK - 12.6941413) * (0.5835786566 - 0.612939346) / (12.84253727 - 12.6941413) 
    
  if( 12.84253727 < T_kK and T_kK < 12.99091615 ) : 
    frac_CaII = 0.5835786566 + (T_kK - 12.84253727) * (0.5536587159 - 0.5835786566) / (12.99091615 - 12.84253727) 
    
  if( 12.99091615 < T_kK and T_kK < 13.12076725 ) : 
    frac_CaII = 0.5536587159 + (T_kK - 12.99091615) * (0.5281195764 - 0.5536587159) / (13.12076725 - 12.99091615) 
    
  if( 13.12076725 < T_kK and T_kK < 13.25059841 ) : 
    frac_CaII = 0.5281195764 + (T_kK - 13.12076725) * (0.5019279773 - 0.5281195764) / (13.25059841 - 13.12076725) 
    
  if( 13.25059841 < T_kK and T_kK < 13.39898584 ) : 
    frac_CaII = 0.5019279773 + (T_kK - 13.25059841) * (0.4722876622 - 0.5019279773) / (13.39898584 - 13.25059841) 
    
  if( 13.39898584 < T_kK and T_kK < 13.54736471 ) : 
    frac_CaII = 0.4722876622 + (T_kK - 13.39898584) * (0.4423677215 - 0.4722876622) / (13.54736471 - 13.39898584) 
    
  if( 13.54736471 < T_kK and T_kK < 13.71448556 ) : 
    frac_CaII = 0.4423677215 + (T_kK - 13.54736471) * (0.4150762615 - 0.4423677215) / (13.71448556 - 13.54736471) 
    
  if( 13.71448556 < T_kK and T_kK < 13.9003159 ) : 
    frac_CaII = 0.4150762615 + (T_kK - 13.71448556) * (0.389350705 - 0.4150762615) / (13.9003159 - 13.71448556) 
    
  if( 13.9003159 < T_kK and T_kK < 14.08615307 ) : 
    frac_CaII = 0.389350705 + (T_kK - 13.9003159) * (0.363848849 - 0.389350705) / (14.08615307 - 13.9003159) 
    
  if( 14.08615307 < T_kK and T_kK < 14.27199025 ) : 
    frac_CaII = 0.363848849 + (T_kK - 14.08615307) * (0.338346993 - 0.363848849) / (14.27199025 - 14.08615307) 
    
  if( 14.27199025 < T_kK and T_kK < 14.45782742 ) : 
    frac_CaII = 0.338346993 + (T_kK - 14.27199025) * (0.312845137 - 0.338346993) / (14.45782742 - 14.27199025) 
    
  if( 14.45782742 < T_kK and T_kK < 14.64356206 ) : 
    frac_CaII = 0.312845137 + (T_kK - 14.45782742) * (0.2839877737 - 0.312845137) / (14.64356206 - 14.45782742) 
    
  if( 14.64356206 < T_kK and T_kK < 14.82928986 ) : 
    frac_CaII = 0.2839877737 + (T_kK - 14.64356206) * (0.2549067098 - 0.2839877737) / (14.82928986 - 14.64356206) 
    
  if( 14.82928986 < T_kK and T_kK < 15.01501767 ) : 
    frac_CaII = 0.2549067098 + (T_kK - 14.82928986) * (0.2258256459 - 0.2549067098) / (15.01501767 - 14.82928986) 
    
  if( 15.01501767 < T_kK and T_kK < 15.20072496 ) : 
    frac_CaII = 0.2258256459 + (T_kK - 15.01501767) * (0.1960734806 - 0.2258256459) / (15.20072496 - 15.01501767) 
    
  if( 15.20072496 < T_kK and T_kK < 15.42373895 ) : 
    frac_CaII = 0.1960734806 + (T_kK - 15.20072496) * (0.1657780426 - 0.1960734806) / (15.42373895 - 15.20072496) 
    
  if( 15.42373895 < T_kK and T_kK < 15.75901514 ) : 
    frac_CaII = 0.1657780426 + (T_kK - 15.42373895) * (0.1450494322 - 0.1657780426) / (15.75901514 - 15.42373895) 
    
  if( 15.75901514 < T_kK and T_kK < 16.16889398 ) : 
    frac_CaII = 0.1450494322 + (T_kK - 15.75901514) * (0.1228827471 - 0.1450494322) / (16.16889398 - 15.75901514) 
    
  if( 16.16889398 < T_kK and T_kK < 16.57895613 ) : 
    frac_CaII = 0.1228827471 + (T_kK - 16.16889398) * (0.1067153025 - 0.1228827471) / (16.57895613 - 16.16889398) 
    
  if( 16.57895613 < T_kK and T_kK < 16.98919229 ) : 
    frac_CaII = 0.1067153025 + (T_kK - 16.57895613) * (0.09624205226 - 0.1067153025) / (16.98919229 - 16.57895613) 
    
  if( 16.98919229 < T_kK and T_kK < 17.39942534 ) : 
    frac_CaII = 0.09624205226 + (T_kK - 16.98919229) * (0.08566711995 - 0.09624205226) / (17.39942534 - 16.98919229) 
    
  if( 17.39942534 < T_kK and T_kK < 17.80959003 ) : 
    frac_CaII = 0.08566711995 + (T_kK - 17.39942534) * (0.07285518272 - 0.08566711995) / (17.80959003 - 17.39942534) 
    
  if( 17.80959003 < T_kK and T_kK < 18.2197423 ) : 
    frac_CaII = 0.07285518272 + (T_kK - 17.80959003) * (0.05963651733 - 0.07285518272) / (18.2197423 - 17.80959003) 
    
  if( 18.2197423 < T_kK and T_kK < 18.62993806 ) : 
    frac_CaII = 0.05963651733 + (T_kK - 18.2197423) * (0.04784140052 - 0.05963651733) / (18.62993806 - 18.2197423) 
    
  if( 18.62993806 < T_kK and T_kK < 19.04029229 ) : 
    frac_CaII = 0.04784140052 + (T_kK - 18.62993806) * (0.04123206782 - 0.04784140052) / (19.04029229 - 18.62993806) 
    
  if( 19.04029229 < T_kK and T_kK < 19.45067138 ) : 
    frac_CaII = 0.04123206782 + (T_kK - 19.04029229) * (0.03543619146 - 0.04123206782) / (19.45067138 - 19.04029229) 
    
  if( 19.45067138 < T_kK and T_kK < 19.86106911 ) : 
    frac_CaII = 0.03543619146 + (T_kK - 19.45067138) * (0.03025040734 - 0.03543619146) / (19.86106911 - 19.45067138) 
    
  if( 19.86106911 < T_kK and T_kK < 20.2715352 ) : 
    frac_CaII = 0.03025040734 + (T_kK - 19.86106911) * (0.02730162814 - 0.03025040734) / (20.2715352 - 19.86106911) 
    
  if( 20.2715352 < T_kK and T_kK < 20.6820075 ) : 
    frac_CaII = 0.02730162814 + (T_kK - 20.2715352) * (0.02455621302 - 0.02730162814) / (20.6820075 - 20.2715352) 
    
  if( 20.6820075 < T_kK and T_kK < 21.09249223 ) : 
    frac_CaII = 0.02455621302 + (T_kK - 20.6820075) * (0.02221752606 - 0.02455621302) / (21.09249223 - 20.6820075) 
    
  if( 21.09249223 < T_kK and T_kK < 21.50297075 ) : 
    frac_CaII = 0.02221752606 + (T_kK - 21.09249223) * (0.01967547503 - 0.02221752606) / (21.50297075 - 21.09249223) 
    
  if( 21.50297075 < T_kK and T_kK < 21.91345548 ) : 
    frac_CaII = 0.01967547503 + (T_kK - 21.50297075) * (0.01733678807 - 0.01967547503) / (21.91345548 - 21.50297075) 
    
  if( 21.91345548 < T_kK and T_kK < 22.32394953 ) : 
    frac_CaII = 0.01733678807 + (T_kK - 21.91345548) * (0.01530314724 - 0.01733678807) / (22.32394953 - 21.91345548) 
    
  if( 22.32394953 < T_kK and T_kK < 22.73445912 ) : 
    frac_CaII = 0.01530314724 + (T_kK - 22.32394953) * (0.01377791662 - 0.01530314724) / (22.73445912 - 22.32394953) 
    
  if( 22.73445912 < T_kK and T_kK < 23.14498735 ) : 
    frac_CaII = 0.01377791662 + (T_kK - 22.73445912) * (0.01286277825 - 0.01377791662) / (23.14498735 - 22.73445912) 
    
  if( 23.14498735 < T_kK and T_kK < 23.55549383 ) : 
    frac_CaII = 0.01286277825 + (T_kK - 23.14498735) * (0.01123586558 - 0.01286277825) / (23.55549383 - 23.14498735) 
    
  if( 23.55549383 < T_kK and T_kK < 23.96603138 ) : 
    frac_CaII = 0.01123586558 + (T_kK - 23.55549383) * (0.01062577333 - 0.01123586558) / (23.96603138 - 23.55549383) 
    
  if( 23.96603138 < T_kK and T_kK < 24.37658136 ) : 
    frac_CaII = 0.01062577333 + (T_kK - 23.96603138) * (0.01042240925 - 0.01062577333) / (24.37658136 - 23.96603138) 
    
  if( 24.37658136 < T_kK and T_kK < 24.7871096 ) : 
    frac_CaII = 0.01042240925 + (T_kK - 24.37658136) * (0.009507270879 - 0.01042240925) / (24.7871096 - 24.37658136) 
    
  if( 24.7871096 < T_kK and T_kK < 25.19763161 ) : 
    frac_CaII = 0.009507270879 + (T_kK - 24.7871096) * (0.008388768422 - 0.009507270879) / (25.19763161 - 24.7871096) 
    
  if( 25.19763161 < T_kK and T_kK < 25.60818781 ) : 
    frac_CaII = 0.008388768422 + (T_kK - 25.19763161) * (0.008388768422 - 0.008388768422) / (25.60818781 - 25.19763161) 
    
  if( 25.60818781 < T_kK and T_kK < 26.01872536 ) : 
    frac_CaII = 0.008388768422 + (T_kK - 25.60818781) * (0.007778676173 - 0.008388768422) / (26.01872536 - 25.60818781) 
    
  if( 26.01872536 < T_kK and T_kK < 26.42926602 ) : 
    frac_CaII = 0.007778676173 + (T_kK - 26.01872536) * (0.007270265966 - 0.007778676173) / (26.42926602 - 26.01872536) 
    
  if( 26.42926602 < T_kK and T_kK < 26.83979736 ) : 
    frac_CaII = 0.007270265966 + (T_kK - 26.42926602) * (0.006456809634 - 0.007270265966) / (26.83979736 - 26.42926602) 
    
  if( 26.83979736 < T_kK and T_kK < 27.25034423 ) : 
    frac_CaII = 0.006456809634 + (T_kK - 26.83979736) * (0.00615176351 - 0.006456809634) / (27.25034423 - 26.83979736) 
    
  if( 27.25034423 < T_kK and T_kK < 27.66090043 ) : 
    frac_CaII = 0.00615176351 + (T_kK - 27.25034423) * (0.00615176351 - 0.00615176351) / (27.66090043 - 27.25034423) 
    
  if( 27.66090043 < T_kK and T_kK < 28.07145662 ) : 
    frac_CaII = 0.00615176351 + (T_kK - 27.66090043) * (0.00615176351 - 0.00615176351) / (28.07145662 - 27.66090043) 
    
  if( 28.07145662 < T_kK and T_kK < 28.48201282 ) : 
    frac_CaII = 0.00615176351 + (T_kK - 28.07145662) * (0.00615176351 - 0.00615176351) / (28.48201282 - 28.07145662) 
    
  if( 28.48201282 < T_kK and T_kK < 28.89255659 ) : 
    frac_CaII = 0.00615176351 + (T_kK - 28.48201282) * (0.005745035344 - 0.00615176351) / (28.89255659 - 28.48201282) 
    
  if( 28.89255659 < T_kK and T_kK < 29.30309103 ) : 
    frac_CaII = 0.005745035344 + (T_kK - 28.89255659) * (0.005033261053 - 0.005745035344) / (29.30309103 - 28.89255659) 
    
  if( 29.30309103 < T_kK and T_kK < 29.71364723 ) : 
    frac_CaII = 0.005033261053 + (T_kK - 29.30309103) * (0.005033261053 - 0.005033261053) / (29.71364723 - 29.30309103) 
    
  if( 29.71364723 < T_kK and T_kK < 29.95624862 ) : 
    frac_CaII = 0.005033261053 + (T_kK - 29.71364723) * (0.005033261053 - 0.005033261053) / (29.95624862 - 29.71364723) 
    
  if( T_kK > 29.95624862 ) : 
    frac_CaII = 0.0 
    
  if( T_kK < 2.326568224 ) : 
    frac_MgII = 1.0 
    
  if( 2.326568224 < T_kK and T_kK < 2.737441152 ) : 
    frac_MgII = 1.046483256 + (T_kK - 2.326568224) * (1.046483256 - 1.046483256) / (2.737441152 - 2.326568224) 
    
  if( 2.737441152 < T_kK and T_kK < 3.14831408 ) : 
    frac_MgII = 1.046483256 + (T_kK - 2.737441152) * (1.046483256 - 1.046483256) / (3.14831408 - 2.737441152) 
    
  if( 3.14831408 < T_kK and T_kK < 3.559187008 ) : 
    frac_MgII = 1.046483256 + (T_kK - 3.14831408) * (1.046483256 - 1.046483256) / (3.559187008 - 3.14831408) 
    
  if( 3.559187008 < T_kK and T_kK < 3.970059936 ) : 
    frac_MgII = 1.046483256 + (T_kK - 3.559187008) * (1.046483256 - 1.046483256) / (3.970059936 - 3.559187008) 
    
  if( 3.970059936 < T_kK and T_kK < 4.380932865 ) : 
    frac_MgII = 1.046483256 + (T_kK - 3.970059936) * (1.046483256 - 1.046483256) / (4.380932865 - 3.970059936) 
    
  if( 4.380932865 < T_kK and T_kK < 4.791805793 ) : 
    frac_MgII = 1.046483256 + (T_kK - 4.380932865) * (1.046483256 - 1.046483256) / (4.791805793 - 4.380932865) 
    
  if( 4.791805793 < T_kK and T_kK < 5.202678721 ) : 
    frac_MgII = 1.046483256 + (T_kK - 4.791805793) * (1.046483256 - 1.046483256) / (5.202678721 - 4.791805793) 
    
  if( 5.202678721 < T_kK and T_kK < 5.613551649 ) : 
    frac_MgII = 1.046483256 + (T_kK - 5.202678721) * (1.046381164 - 1.046483256) / (5.613551649 - 5.202678721) 
    
  if( 5.613551649 < T_kK and T_kK < 6.024424577 ) : 
    frac_MgII = 1.046381164 + (T_kK - 5.613551649) * (1.045258151 - 1.046381164) / (6.024424577 - 5.613551649) 
    
  if( 6.024424577 < T_kK and T_kK < 6.435297505 ) : 
    frac_MgII = 1.045258151 + (T_kK - 6.024424577) * (1.04372677 - 1.045258151) / (6.435297505 - 6.024424577) 
    
  if( 6.435297505 < T_kK and T_kK < 6.846170434 ) : 
    frac_MgII = 1.04372677 + (T_kK - 6.435297505) * (1.041378653 - 1.04372677) / (6.846170434 - 6.435297505) 
    
  if( 6.846170434 < T_kK and T_kK < 7.257043362 ) : 
    frac_MgII = 1.041378653 + (T_kK - 6.846170434) * (1.03974518 - 1.041378653) / (7.257043362 - 6.846170434) 
    
  if( 7.257043362 < T_kK and T_kK < 7.66791629 ) : 
    frac_MgII = 1.03974518 + (T_kK - 7.257043362) * (1.038315891 - 1.03974518) / (7.66791629 - 7.257043362) 
    
  if( 7.66791629 < T_kK and T_kK < 8.078789218 ) : 
    frac_MgII = 1.038315891 + (T_kK - 7.66791629) * (1.035048945 - 1.038315891) / (8.078789218 - 7.66791629) 
    
  if( 8.078789218 < T_kK and T_kK < 8.489662146 ) : 
    frac_MgII = 1.035048945 + (T_kK - 8.078789218) * (1.031271538 - 1.035048945) / (8.489662146 - 8.078789218) 
    
  if( 8.489662146 < T_kK and T_kK < 8.900535075 ) : 
    frac_MgII = 1.031271538 + (T_kK - 8.489662146) * (1.026881579 - 1.031271538) / (8.900535075 - 8.489662146) 
    
  if( 8.900535075 < T_kK and T_kK < 9.311408003 ) : 
    frac_MgII = 1.026881579 + (T_kK - 8.900535075) * (1.021164424 - 1.026881579) / (9.311408003 - 8.900535075) 
    
  if( 9.311408003 < T_kK and T_kK < 9.722280931 ) : 
    frac_MgII = 1.021164424 + (T_kK - 9.311408003) * (1.011057309 - 1.021164424) / (9.722280931 - 9.311408003) 
    
  if( 9.722280931 < T_kK and T_kK < 10.13315386 ) : 
    frac_MgII = 1.011057309 + (T_kK - 9.722280931) * (1.000439734 - 1.011057309) / (10.13315386 - 9.722280931) 
    
  if( 10.13315386 < T_kK and T_kK < 10.54402679 ) : 
    frac_MgII = 1.000439734 + (T_kK - 10.13315386) * (0.9902305272 - 1.000439734) / (10.54402679 - 10.13315386) 
    
  if( 10.54402679 < T_kK and T_kK < 10.95489972 ) : 
    frac_MgII = 0.9902305272 + (T_kK - 10.54402679) * (0.9752229933 - 0.9902305272) / (10.95489972 - 10.54402679) 
    
  if( 10.95489972 < T_kK and T_kK < 11.29106847 ) : 
    frac_MgII = 0.9752229933 + (T_kK - 10.95489972) * (0.9585674017 - 0.9752229933) / (11.29106847 - 10.95489972) 
    
  if( 11.29106847 < T_kK and T_kK < 11.53385702 ) : 
    frac_MgII = 0.9585674017 + (T_kK - 11.29106847) * (0.9311872815 - 0.9585674017) / (11.53385702 - 11.29106847) 
    
  if( 11.53385702 < T_kK and T_kK < 11.75796953 ) : 
    frac_MgII = 0.9311872815 + (T_kK - 11.53385702) * (0.9038606382 - 0.9311872815) / (11.75796953 - 11.53385702) 
    
  if( 11.75796953 < T_kK and T_kK < 11.96340599 ) : 
    frac_MgII = 0.9038606382 + (T_kK - 11.75796953) * (0.8789297553 - 0.9038606382) / (11.96340599 - 11.75796953) 
    
  if( 11.96340599 < T_kK and T_kK < 12.13149037 ) : 
    frac_MgII = 0.8789297553 + (T_kK - 11.96340599) * (0.8519212989 - 0.8789297553) / (12.13149037 - 11.96340599) 
    
  if( 12.13149037 < T_kK and T_kK < 12.28089871 ) : 
    frac_MgII = 0.8519212989 + (T_kK - 12.13149037) * (0.8244074868 - 0.8519212989) / (12.28089871 - 12.13149037) 
    
  if( 12.28089871 < T_kK and T_kK < 12.43030705 ) : 
    frac_MgII = 0.8244074868 + (T_kK - 12.28089871) * (0.7985781937 - 0.8244074868) / (12.43030705 - 12.28089871) 
    
  if( 12.43030705 < T_kK and T_kK < 12.57971539 ) : 
    frac_MgII = 0.7985781937 + (T_kK - 12.43030705) * (0.7716258879 - 0.7985781937) / (12.57971539 - 12.43030705) 
    
  if( 12.57971539 < T_kK and T_kK < 12.72912372 ) : 
    frac_MgII = 0.7716258879 + (T_kK - 12.57971539) * (0.7446735821 - 0.7716258879) / (12.72912372 - 12.57971539) 
    
  if( 12.72912372 < T_kK and T_kK < 12.85985602 ) : 
    frac_MgII = 0.7446735821 + (T_kK - 12.72912372) * (0.7193122111 - 0.7446735821) / (12.85985602 - 12.72912372) 
    
  if( 12.85985602 < T_kK and T_kK < 12.99058831 ) : 
    frac_MgII = 0.7193122111 + (T_kK - 12.85985602) * (0.6921727365 - 0.7193122111) / (12.99058831 - 12.85985602) 
    
  if( 12.99058831 < T_kK and T_kK < 13.12132061 ) : 
    frac_MgII = 0.6921727365 + (T_kK - 12.99058831) * (0.6642845868 - 0.6921727365) / (13.12132061 - 12.99058831) 
    
  if( 13.12132061 < T_kK and T_kK < 13.2520529 ) : 
    frac_MgII = 0.6642845868 + (T_kK - 13.12132061) * (0.6374258654 - 0.6642845868) / (13.2520529 - 13.12132061) 
    
  if( 13.2520529 < T_kK and T_kK < 13.3827852 ) : 
    frac_MgII = 0.6374258654 + (T_kK - 13.2520529) * (0.6096313 - 0.6374258654) / (13.3827852 - 13.2520529) 
    
  if( 13.3827852 < T_kK and T_kK < 13.49484145 ) : 
    frac_MgII = 0.6096313 + (T_kK - 13.3827852) * (0.5864223701 - 0.6096313) / (13.49484145 - 13.3827852) 
    
  if( 13.49484145 < T_kK and T_kK < 13.60689771 ) : 
    frac_MgII = 0.5864223701 + (T_kK - 13.49484145) * (0.5624647649 - 0.5864223701) / (13.60689771 - 13.49484145) 
    
  if( 13.60689771 < T_kK and T_kK < 13.71895396 ) : 
    frac_MgII = 0.5624647649 + (T_kK - 13.60689771) * (0.5381328222 - 0.5624647649) / (13.71895396 - 13.60689771) 
    
  if( 13.71895396 < T_kK and T_kK < 13.83101021 ) : 
    frac_MgII = 0.5381328222 + (T_kK - 13.71895396) * (0.5130522043 - 0.5381328222) / (13.83101021 - 13.71895396) 
    
  if( 13.83101021 < T_kK and T_kK < 13.94306647 ) : 
    frac_MgII = 0.5130522043 + (T_kK - 13.83101021) * (0.4887202616 - 0.5130522043) / (13.94306647 - 13.83101021) 
    
  if( 13.94306647 < T_kK and T_kK < 14.05512272 ) : 
    frac_MgII = 0.4887202616 + (T_kK - 13.94306647) * (0.4636396437 - 0.4887202616) / (14.05512272 - 13.94306647) 
    
  if( 14.05512272 < T_kK and T_kK < 14.16717897 ) : 
    frac_MgII = 0.4636396437 + (T_kK - 14.05512272) * (0.4381846882 - 0.4636396437) / (14.16717897 - 14.05512272) 
    
  if( 14.16717897 < T_kK and T_kK < 14.27923523 ) : 
    frac_MgII = 0.4381846882 + (T_kK - 14.16717897) * (0.4134784079 - 0.4381846882) / (14.27923523 - 14.16717897) 
    
  if( 14.27923523 < T_kK and T_kK < 14.39129148 ) : 
    frac_MgII = 0.4134784079 + (T_kK - 14.27923523) * (0.3876491149 - 0.4134784079) / (14.39129148 - 14.27923523) 
    
  if( 14.39129148 < T_kK and T_kK < 14.52202377 ) : 
    frac_MgII = 0.3876491149 + (T_kK - 14.39129148) * (0.3614454843 - 0.3876491149) / (14.52202377 - 14.39129148) 
    
  if( 14.52202377 < T_kK and T_kK < 14.67143211 ) : 
    frac_MgII = 0.3614454843 + (T_kK - 14.52202377) * (0.3333701657 - 0.3614454843) / (14.67143211 - 14.52202377) 
    
  if( 14.67143211 < T_kK and T_kK < 14.82084045 ) : 
    frac_MgII = 0.3333701657 + (T_kK - 14.67143211) * (0.30641786 - 0.3333701657) / (14.82084045 - 14.67143211) 
    
  if( 14.82084045 < T_kK and T_kK < 14.97024879 ) : 
    frac_MgII = 0.30641786 + (T_kK - 14.82084045) * (0.2794655542 - 0.30641786) / (14.97024879 - 14.82084045) 
    
  if( 14.97024879 < T_kK and T_kK < 15.11965712 ) : 
    frac_MgII = 0.2794655542 + (T_kK - 14.97024879) * (0.2525132484 - 0.2794655542) / (15.11965712 - 14.97024879) 
    
  if( 15.11965712 < T_kK and T_kK < 15.26906546 ) : 
    frac_MgII = 0.2525132484 + (T_kK - 15.11965712) * (0.2255609426 - 0.2525132484) / (15.26906546 - 15.11965712) 
    
  if( 15.26906546 < T_kK and T_kK < 15.45582588 ) : 
    frac_MgII = 0.2255609426 + (T_kK - 15.26906546) * (0.1976727929 - 0.2255609426) / (15.45582588 - 15.26906546) 
    
  if( 15.45582588 < T_kK and T_kK < 15.69861443 ) : 
    frac_MgII = 0.1976727929 + (T_kK - 15.45582588) * (0.1734210653 - 0.1976727929) / (15.69861443 - 15.45582588) 
    
  if( 15.69861443 < T_kK and T_kK < 15.96007902 ) : 
    frac_MgII = 0.1734210653 + (T_kK - 15.69861443) * (0.1474313419 - 0.1734210653) / (15.96007902 - 15.69861443) 
    
  if( 15.96007902 < T_kK and T_kK < 16.22154361 ) : 
    frac_MgII = 0.1474313419 + (T_kK - 15.96007902) * (0.1214416185 - 0.1474313419) / (16.22154361 - 15.96007902) 
    
  if( 16.22154361 < T_kK and T_kK < 16.55771237 ) : 
    frac_MgII = 0.1214416185 + (T_kK - 16.22154361) * (0.0970270298 - 0.1214416185) / (16.55771237 - 16.22154361) 
    
  if( 16.55771237 < T_kK and T_kK < 16.9685853 ) : 
    frac_MgII = 0.0970270298 + (T_kK - 16.55771237) * (0.07681280046 - 0.0970270298) / (16.9685853 - 16.55771237) 
    
  if( 16.9685853 < T_kK and T_kK < 17.37945823 ) : 
    frac_MgII = 0.07681280046 + (T_kK - 16.9685853) * (0.05731321559 - 0.07681280046) / (17.37945823 - 16.9685853) 
    
  if( 17.37945823 < T_kK and T_kK < 17.79033116 ) : 
    frac_MgII = 0.05731321559 + (T_kK - 17.37945823) * (0.04567471991 - 0.05731321559) / (17.79033116 - 17.37945823) 
    
  if( 17.79033116 < T_kK and T_kK < 18.20120409 ) : 
    frac_MgII = 0.04567471991 + (T_kK - 17.79033116) * (0.03628224972 - 0.04567471991) / (18.20120409 - 17.79033116) 
    
  if( 18.20120409 < T_kK and T_kK < 18.61207701 ) : 
    frac_MgII = 0.03628224972 + (T_kK - 18.20120409) * (0.02770651606 - 0.03628224972) / (18.61207701 - 18.20120409) 
    
  if( 18.61207701 < T_kK and T_kK < 19.02294994 ) : 
    frac_MgII = 0.02770651606 + (T_kK - 18.61207701) * (0.02290818889 - 0.02770651606) / (19.02294994 - 18.61207701) 
    
  if( 19.02294994 < T_kK and T_kK < 19.43382287 ) : 
    frac_MgII = 0.02290818889 + (T_kK - 19.02294994) * (0.01851822999 - 0.02290818889) / (19.43382287 - 19.02294994) 
    
  if( 19.43382287 < T_kK and T_kK < 19.8446958 ) : 
    frac_MgII = 0.01851822999 + (T_kK - 19.43382287) * (0.01453663937 - 0.01851822999) / (19.8446958 - 19.43382287) 
    
  if( 19.8446958 < T_kK and T_kK < 20.25556873 ) : 
    frac_MgII = 0.01453663937 + (T_kK - 19.8446958) * (0.01229061388 - 0.01453663937) / (20.25556873 - 19.8446958) 
    
  if( 20.25556873 < T_kK and T_kK < 20.66644165 ) : 
    frac_MgII = 0.01229061388 + (T_kK - 20.25556873) * (0.01024877254 - 0.01229061388) / (20.66644165 - 20.25556873) 
    
  if( 20.66644165 < T_kK and T_kK < 21.07731458 ) : 
    frac_MgII = 0.01024877254 + (T_kK - 20.66644165) * (0.008104839123 - 0.01024877254) / (21.07731458 - 20.66644165) 
    
  if( 21.07731458 < T_kK and T_kK < 21.48818751 ) : 
    frac_MgII = 0.008104839123 + (T_kK - 21.07731458) * (0.007288102584 - 0.008104839123) / (21.48818751 - 21.07731458) 
    
  if( 21.48818751 < T_kK and T_kK < 21.89906044 ) : 
    frac_MgII = 0.007288102584 + (T_kK - 21.48818751) * (0.005756721574 - 0.007288102584) / (21.89906044 - 21.48818751) 
    
  if( 21.89906044 < T_kK and T_kK < 22.30993337 ) : 
    frac_MgII = 0.005756721574 + (T_kK - 21.89906044) * (0.005042077102 - 0.005756721574) / (22.30993337 - 21.89906044) 
    
  if( 22.30993337 < T_kK and T_kK < 22.72080629 ) : 
    frac_MgII = 0.005042077102 + (T_kK - 22.30993337) * (0.004327432631 - 0.005042077102) / (22.72080629 - 22.30993337) 
    
  if( 22.72080629 < T_kK and T_kK < 23.13167922 ) : 
    frac_MgII = 0.004327432631 + (T_kK - 22.72080629) * (0.003408604024 - 0.004327432631) / (23.13167922 - 22.72080629) 
    
  if( 23.13167922 < T_kK and T_kK < 23.54255215 ) : 
    frac_MgII = 0.003408604024 + (T_kK - 23.13167922) * (0.00320441989 - 0.003408604024) / (23.54255215 - 23.13167922) 
    
  if( 23.54255215 < T_kK and T_kK < 23.95342508 ) : 
    frac_MgII = 0.00320441989 + (T_kK - 23.54255215) * (0.00279605162 - 0.00320441989) / (23.95342508 - 23.54255215) 
    
  if( 23.95342508 < T_kK and T_kK < 24.36429801 ) : 
    frac_MgII = 0.00279605162 + (T_kK - 23.95342508) * (0.003000235755 - 0.00279605162) / (24.36429801 - 23.95342508) 
    
  if( 24.36429801 < T_kK and T_kK < 24.77517094 ) : 
    frac_MgII = 0.003000235755 + (T_kK - 24.36429801) * (0.002081407148 - 0.003000235755) / (24.77517094 - 24.36429801) 
    
  if( 24.77517094 < T_kK and T_kK < 25.18604386 ) : 
    frac_MgII = 0.002081407148 + (T_kK - 24.77517094) * (0.002081407148 - 0.002081407148) / (25.18604386 - 24.77517094) 
    
  if( 25.18604386 < T_kK and T_kK < 25.59691679 ) : 
    frac_MgII = 0.002081407148 + (T_kK - 25.18604386) * (0.001877223014 - 0.002081407148) / (25.59691679 - 25.18604386) 
    
  if( 25.59691679 < T_kK and T_kK < 26.00778972 ) : 
    frac_MgII = 0.001877223014 + (T_kK - 25.59691679) * (0.00126467061 - 0.001877223014) / (26.00778972 - 25.59691679) 
    
  if( 26.00778972 < T_kK and T_kK < 26.41866265 ) : 
    frac_MgII = 0.00126467061 + (T_kK - 26.00778972) * (0.001162578542 - 0.00126467061) / (26.41866265 - 26.00778972) 
    
  if( 26.41866265 < T_kK and T_kK < 26.82953558 ) : 
    frac_MgII = 0.001162578542 + (T_kK - 26.41866265) * (0.001060486475 - 0.001162578542) / (26.82953558 - 26.41866265) 
    
  if( 26.82953558 < T_kK and T_kK < 27.2404085 ) : 
    frac_MgII = 0.001060486475 + (T_kK - 26.82953558) * (0.0009583944075 - 0.001060486475) / (27.2404085 - 26.82953558) 
    
  if( 27.2404085 < T_kK and T_kK < 27.65128143 ) : 
    frac_MgII = 0.0009583944075 + (T_kK - 27.2404085) * (0.0009583944075 - 0.0009583944075) / (27.65128143 - 27.2404085) 
    
  if( 27.65128143 < T_kK and T_kK < 28.06215436 ) : 
    frac_MgII = 0.0009583944075 + (T_kK - 27.65128143) * (0.0009583944075 - 0.0009583944075) / (28.06215436 - 27.65128143) 
    
  if( 28.06215436 < T_kK and T_kK < 28.47302729 ) : 
    frac_MgII = 0.0009583944075 + (T_kK - 28.06215436) * (0.0009583944075 - 0.0009583944075) / (28.47302729 - 28.06215436) 
    
  if( 28.47302729 < T_kK and T_kK < 28.88390022 ) : 
    frac_MgII = 0.0009583944075 + (T_kK - 28.47302729) * (0.0009583944075 - 0.0009583944075) / (28.88390022 - 28.47302729) 
    
  if( 28.88390022 < T_kK and T_kK < 29.29477315 ) : 
    frac_MgII = 0.0009583944075 + (T_kK - 28.88390022) * (-6.252626616e-05 - 0.0009583944075) / (29.29477315 - 28.88390022) 
    
  if( 29.29477315 < T_kK and T_kK < 29.70564607 ) : 
    frac_MgII = -6.252626616e-05 + (T_kK - 29.29477315) * (-0.0001646183335 - -6.252626616e-05) / (29.70564607 - 29.29477315) 
    
  if( 29.70564607 < T_kK and T_kK < 29.94843462 ) : 
    frac_MgII = -0.0001646183335 + (T_kK - 29.70564607) * (0.0009583944075 - -0.0001646183335) / (29.94843462 - 29.70564607) 
    
  if( T_kK > 29.94843462 ) : 
    frac_MgII = 0.0 
    
  Q_H = (10.0**log10Q_HI_Lyman) + (10.0**log10Q_HI_neutral)
  Q_Mg = (10.0**log10Q_MgII)
  Q_Ca = (10.0**log10Q_CaII)

  if T_kK < 4.:
	Q_H = 0.0

  if T_kK < 5.:
	Q_Mg = 0.0


  if esc_CaII > 1.0:
	esc_CaII = 1.0
  if esc_MgII > 1.0:
	esc_MgII = 1.0
  if esc_HI > 1.0:
	esc_HI = 1.0

  if frac_HI > 1.0:
	frac_HI = 1.0
  if frac_MgII > 1.0:
	frac_MgII = 1.0
  if frac_CaII > 1.0:
	frac_CaII = 1.0

  return Q_H, Q_Mg, Q_Ca, esc_HI, esc_MgII, esc_CaII, frac_HI, frac_MgII, frac_CaII



def transform_SI_to_CGS(mesh_arr, data_arr, ws, hs, nnodes, coordinates, connectivity):
	Bx_arr = data_arr[0]
        By_arr = data_arr[1]
        Bz_arr = data_arr[2]
        Ex_arr = data_arr[3]
        Ey_arr = data_arr[4]
        Ez_arr = data_arr[5]
        psi_arr = data_arr[6]
        phi_arr = data_arr[7]
        rhoi_arr = data_arr[8]
        rhon_arr = data_arr[9]
        Vxi_arr = data_arr[10]
        Vyi_arr = data_arr[11]
        Vxn_arr = data_arr[12]
        Vyn_arr = data_arr[13]
        Ti_arr = data_arr[14]
        Tn_arr = data_arr[15]	

	for i in range(0, len(data_arr[0])):
		#Rhos
		data_arr[8][i] = data_arr[8][i]*0.001
		data_arr[9][i] = data_arr[9][i]*0.001

		#Velocities
		data_arr[10][i] = data_arr[10][i]*100.0
		data_arr[11][i] = data_arr[11][i]*100.0
		data_arr[12][i] = data_arr[12][i]*100.0
		data_arr[13][i] = data_arr[13][i]*100.0

	#Mesh and coordinates
	for i in range(0, len(mesh_arr)):
		mesh_arr[i][0] = mesh_arr[i][0]*100.
		mesh_arr[i][1] = mesh_arr[i][1]*100.

	
	for i in range(0, len(coordinates)):
		coordinates[i][0] = coordinates[i][0]*100.
		coordinates[i][1] = coordinates[i][1]*100.

	for i in range(0, len(hs)):
		hs[i] = 100.0*hs[i]

	for i in range(0, len(ws)):
		ws[i] = 100.0*ws[i]

	return mesh_arr, data_arr, ws, hs, nnodes, coordinates, connectivity



#This is the input file:

mhdin_filename = str((sys.argv)[1])
plt_filename = str((sys.argv)[2])
MHDinfile = open(mhdin_filename, "r")

imax = 149
jmax = 299


#Get CFmesh data, without any unit transformations (units in SI):
mesh_arr, data_arr, ws, hs, nnodes, coordinates, connectivity = read_CFmesh_data(MHDinfile)

#Transform to SI:
mesh_arr, data_arr, ws, hs, nnodes, coordinates, connectivity = transform_SI_to_CGS(mesh_arr, data_arr, ws, hs, nnodes, coordinates, connectivity)


#Sort cells:
cells_arr = sort_cells(imax, jmax)



#Compute radiative flux:
Eloss_arr, Qlocal, Qcoronal, Qphotospheric, taus, mcs,  Qrad_H, Qrad_Mg, Qrad_Ca = computeRadiation(mesh_arr, cells_arr, data_arr, ws, hs, imax, jmax)



Ts = data_arr[14]

rhos = []
for i in range(0, len(data_arr[9])):
	rhos.append(data_arr[9][i] + data_arr[8][i])


filename = plt_filename 

writePLT(mesh_arr[:][0], mesh_arr[:][1], Qlocal, Qcoronal, Qphotospheric, Ts, rhos, taus, mcs,  Qrad_H, Qrad_Mg, Qrad_Ca, filename, nnodes, coordinates, connectivity)


#for j in range(0, imax*jmax):
#	if cells_arr[j][0] == 50:
#		print mesh_arr[j][1], Eloss_arr[j]




