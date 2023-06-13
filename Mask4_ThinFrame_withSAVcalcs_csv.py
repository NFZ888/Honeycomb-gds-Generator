import numpy as np
import gdspy as gd
import math
import csv

#PyScript to automatically generate parameter sweep .gds files (used for LMIS1_V4 Mask)
#Created by Nicolas Zaugg as part of Semester Project at LMIS1, EPFL
#Uses gdspy: https://gdspy.readthedocs.io/en/stable/gettingstarted.html#loading-a-gdsii-file

#Parameters of the unit honecomb cell
theta = 30
l_b = 60
t = 10
fillet_r = 8
t_inc = 5

#Parameters of the honeycomb array
device_size = 4000

#paramter sweep parameters
psa_s = 5500
global_offset_y = 0



# The GDSII file is called a library, which contains multiple cells.
lib = gd.GdsLibrary()

#full_gds = lib.new_cell("Full_Gds")
device_array = lib.new_cell("Device_Array")

#3
for j in range(3):

    #6
    for i in range(7):

        #find position/offset for parameter sweep array
        sweep_offset_y = -2*(math.floor(i/10)*psa_s+1)
        sweep_offset_x = int(str(i)[-1])*psa_s

        if (i>2): 
            offset_x = offset_x - 3*psa_s

        #caluclate outer length
        l = l_b + (t/2) * (1/math.cos(math.radians(theta)))

        #trig
        c = math.cos(math.radians(60))
        s = math.sin(math.radians(60))
        offset_x = t/2
        offset_y = (l-l_b)

        #Width Unit cell
        w = t+2*(math.sin(math.radians(60))*l_b)
        #Height Unit cell
        h = 2*c*l+l

        #inner wdith, height
        w_b = 2*(math.sin(math.radians(60))*l_b)
        h_b = 2*c*l_b+l_b

        #nb of cells
        n_rows = math.ceil(0.5*math.ceil(device_size/h))+1
        n_columns = math.floor(device_size/w)


        #Unit Hex
        hex = lib.new_cell('Hex'+str(j)+str(i))

        #Inner & outer hex geometry parametrization
        points_ohex = [(0, l*c), (0, l*c+l), (l*s, 2*l*c+l), (2*l*s, l*c+l), (2*l*s, l*c), (l*s, 0)]
        points_ihex = [(0, l_b*c), (0, l_b*c+l_b), (l_b*s, 2*l_b*c+l_b), (2*l_b*s, l_b*c+l_b), (2*l_b*s, l_b*c), (l_b*s, 0)]
        hex_2D = gd.boolean(gd.Polygon(points_ohex), gd.Polygon(points_ihex).translate(offset_x, offset_y).fillet(fillet_r), "not")

        hex.add(hex_2D)

        #Array geometry parametrization
        hex_array = lib.new_cell("Hex_Array"+str(j)+str(i))
        
        #Offset rows
        magic_offset_x = 500
        offset_rows = (w/2-magic_offset_x, -(h/2 + l/2))

        #Arrays of even and odd rows
        arr1 = gd.CellArray(hex, n_columns, n_rows, (w, h+l), (-magic_offset_x,0))
        arr2 = gd.CellArray(hex, n_columns, n_rows, (w, h+l), (offset_rows))

        #frame
        iframe = gd.Rectangle([t, t], [2500-t, 2500-t])
        oframe = gd.Rectangle([0, 0], [2500, 2500] )
        frame = gd.boolean(oframe, iframe, 'not')

        #cut
        arr1_Cut = gd.boolean(arr1, iframe, 'and')
        arr2_Cut = gd.boolean(arr2, iframe, 'and')
        hex_array.add(arr1_Cut)
        hex_array.add(arr2_Cut)
        hex_array.add(frame)

        #refsquares
        if (i == 6 ):
            ref_sq = lib.new_cell('Ref_sq'+str(j)+str(i))
            ref_sq.add(oframe)
            device_array.add(gd.CellArray(ref_sq, 1, 2, (psa_s, psa_s), (sweep_offset_x, sweep_offset_y+global_offset_y)))
        else:
            device_array.add(gd.CellArray(hex_array, 1, 2, (psa_s, psa_s), (sweep_offset_x, sweep_offset_y+global_offset_y)))

        #SA/V caculations, outputs to a csv file 
        #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        n_x =(2500/w)
        n_y = (2500/(h_b+l_b+2*(math.sqrt(t*t-(t/2)*(t/2)))))*2
        area = hex_array.area()
        inner_path = n_x*n_y*6*l_b
        with open("SAV.csv", 'a', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow([l_b, t, ((inner_path*100+4*2500*100+2*area)/area*100)])
        #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        if(i > 1):
            t = t*2
        #new line thickness
        else:
            t = t + t_inc

    #Housekeeping for next loop
    t=10
    global_offset_y = global_offset_y + sweep_offset_y - 10918
    l_b = l_b + 20
    fillet_r = fillet_r + 3

lib.write_gds("TEST.gds")