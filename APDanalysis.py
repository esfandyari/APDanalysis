import sys
import openpyxl as oxl
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import peakutils
from scipy import interpolate
from bpeak_v2 import *

#################################### parameter of analysis ##############################
minimum_peak= .02		        ###minimum height of peaks
delta   = .003			        #root finding initial interval
tol     = 0.000001                      #tolerance of the root finding function
max_it  = 500                           ###maximum iteration to find a root>>> increase to enhance the accuracy of APDs, would take a bit longer
b_tol   = 0.0000000001                  #tolerance of base-line finding function
p_enhance 	= 5			###enhancemet parameter for peak enhancer>>> compares this many data point to find the peak
time_lapse 	= .08                   ##time steps of raw data>>> time lapse, should be adjusted for FRET or ANEPPS
##################################### main ##############################################

#assign the execution arguments
sys_arg     = sys.argv

#import data from excel file
data        = importer(str(sys_arg[1])+'.xlsx','','A','B') 							
print 'file is loaded'

#finds the baseline of the raw data
baseline    = peakutils.baseline(data[1][:], deg=2, max_it=10000, tol=b_tol)

plot1 = plt.plot(data[0][:],baseline)
plot2 = plt.plot(data[0][:],data[1][:])
plt.show()

# the interpolationg function of baseline
base_f      = interpolate.interp1d(data[0][:],baseline)

# subtracting baseline from data
data        = data_subtract_baseline(data,base_f)

#calculates the time step
time_laps  = data[0][8] - data[0][7]

print("The time steps are %f\n" % time_laps)

#regenerate the baseline and its interpolating function
baseline    = peakutils.baseline(data[1][:], deg=1, max_it=10000, tol=b_tol)
base_f      = interpolate.interp1d(data[0][:],baseline)


data_f      = interpolate.interp1d(data[0][:],data[1][:]) 						# return data as interpolationg function
###################################################################################

peak_id     = signal.find_peaks_cwt(data[1][:], np.arange(1,10))                                        #finding the index of the peaks
peak_h      = peak_height(peak_id,base_f,data) 
peak_id     = peak_filter(peak_h,peak_id,minimum_peak*.8)
peak_id     = peak_modifier(peak_id,data,p_enhance) 							#enhance the peaks
peak_id     = np.unique(peak_id)			
peak_h      = peak_height(peak_id,base_f,data) 								#return the height of peaks
peak_id     = peak_filter(peak_h,peak_id,minimum_peak)							#eliminate the peaks shorter than minimum peaks
peak_h      = peak_height(peak_id,base_f,data)

plot1 = plt.plot(data[0][:],data[1][:])
plot2 = plt.plot(data[0][:],baseline)
#plot10= plt.plot(data[0][:]+.5,data_f(data[0][:])+.01,'--b')
plot3 = plt.plot(data[0][peak_id],data[1][peak_id],'ro')
plt.show()
							#return the height after elimination
RR              = rr_calc(data,peak_id)
#return the RR
'''
#peak_s_p        = peak_split_point(data,peak_id,RR) 							#return the splitting points
#peak_data       = peak_data_gen(peak_s_p,data,data_f,time_lapse) 					#return a matrix of each data- a splitted peck on each row
#peak_location   = (peak_locator(data)).astype(int) 						#return the time of each peak for each splitted peak
#peak_data_av    = data_averager(peak_data,peak_location) 						# work with one peak average
#data[1][:]      = peak_data_av
'''
######################################################################################

figure1     = plt.plot(data[0][:],data[1][:],'.-b')
plt.tight_layout()
plt.show()
answer 		= 0
var 		= 0
while(answer) !=1:
        plt.clf()
	var             = float(raw_input("Please enter baseline y: "))
	baseline        = np.full(len(data[1]), var, dtype=np.float)

	# redefines both baseline and data points by interpolating function
        data_f          = interpolate.interp1d(data[0][:],data[1][:])
        base_f          = interpolate.interp1d(data[0][:],baseline)
        """
        # locates the peaks
        peak_id_temp    = signal.find_peaks_cwt(data[1][:], np.arange(5,10))
        
        # modifies the peak by vicinity search
        peak_id         = peak_modifier(peak_id_temp,data,7)
        
        # filters the peaks with respect to the minimum peak height
        peak_h_temp     = peak_height(peak_id,base_f,data)
        peak_id         = peak_filter(peak_h_temp,peak_id,minimum_peak)
        """
        peak_h          = peak_height(peak_id,base_f,data)
        
        # finds the 20,50,70 and 90% of the peaks height
        p_h_90          = peak_h_90(peak_h)
        p_h_70          = peak_h_70(peak_h)
        p_h_50          = peak_h_50(peak_h)
        p_h_20          = peak_h_20(peak_h)

        # cross-line of the peaks at 50% and 90%
        p_x_90rl        = peak_r_l(data,baseline,p_h_90,peak_id,delta,tol,data_f,max_it)
        p_x_70rl        = peak_r_l(data,baseline,p_h_70,peak_id,delta,tol,data_f,max_it)
        p_x_50rl        = peak_r_l(data,baseline,p_h_50,peak_id,delta,tol,data_f,max_it)
        p_x_20rl        = peak_r_l(data,baseline,p_h_20,peak_id,delta,tol,data_f,max_it)

        # plots
        plot1 = plt.plot(data[0][:],data[1][:])
        plot2 = plt.plot(data[0][:],baseline)
        #plot10= plt.plot(data[0][:]+.5,data_f(data[0][:])+.01,'--b')
        plot3 = plt.plot(data[0][peak_id],data[1][peak_id],'ro')

        plot5 = plt.plot(p_x_90rl[0], data_f(p_x_90rl[0]),'go')
        plot6 = plt.plot(p_x_70rl[0], data_f(p_x_70rl[0]),'bo')
        plot7 = plt.plot(p_x_50rl[0], data_f(p_x_50rl[0]),'go')
        plot8 = plt.plot(p_x_20rl[0], data_f(p_x_20rl[0]),'yo')
        plot9 = plt.plot(p_x_50rl[1], data_f(p_x_50rl[0]),'go')
        plot10 = plt.plot(p_x_90rl[1], data_f(p_x_90rl[0]),'go')
        plot11 = plt.plot(p_x_70rl[1], data_f(p_x_70rl[0]),'bo')
        plot12 = plt.plot(p_x_20rl[1], data_f(p_x_20rl[0]),'yo')
        figure3     = plt.plot(data[0][:],data[1][:],'.-b')
        plt.tight_layout()
        plt.show(block=False)
        answer      = int(raw_input("1 to continue 0 to repeat: "))

#######################################################################################
plt.clf()
figure1     = plt.plot(data[0][:],data[1][:],'.-b')
plt.tight_layout()
#plt.show()
# calculates the properties
QT_90           = np.abs(p_x_90rl[0] -  p_x_90rl[1])
QT_70           = np.abs(p_x_70rl[0] -  p_x_70rl[1])
QT_50           = np.abs(p_x_50rl[0] -  p_x_50rl[1])
QT_20           = np.abs(p_x_20rl[0] -  p_x_20rl[1])
RR_av           = np.sum(RR)/len(RR)
RR_av_Sqrt      = np.sqrt(RR_av)
av_QT90         = np.sum(QT_90)/len(QT_90)
av_QT70         = np.sum(QT_70)/len(QT_70)
av_QT50         = np.sum(QT_50)/len(QT_50)
av_QT20         = np.sum(QT_20)/len(QT_20)
QTc_90          = (av_QT90) / RR_av_Sqrt
QTc_70          = (av_QT70) / RR_av_Sqrt
QTc_50          = (av_QT50) / RR_av_Sqrt
QTc_20          = (av_QT20) / RR_av_Sqrt
BPM             = 60 * len(RR)/ np.sum(RR)

################################# slope of gradual depolarization ################
'''
answer = 0
x1_in     = intersect_left(data[0][peak_id[0]],base_f,data_f,delta,tol)
y1_in     = var
while(answer) !=1:
        var             = float(raw_input("Please enter 2nd baseline y: "))
        baseline_2nd    = np.full(len(data[1]), var, dtype=np.float)
        base_f_2nd      = interpolate.interp1d(data[0][:],baseline_2nd)
        x2_in           = intersect_left(data[0][peak_id[0]],base_f_2nd,data_f,delta,tol)
        y2_in           = var
        print(str(y2_in),'\t' + str(y1_in),'\t' +str(x2_in),'\t' +str(x1_in))
        slope           =(y2_in - y1_in)/(x2_in - x1_in)
        # plots
        plot1 = plt.plot(data[0][:],data[1][:])
        plot2 = plt.plot(data[0][:],baseline_2nd)
        plot3 = plt.plot(data[0][peak_id],data[1][peak_id],'ro')
        figure3     = plt.plot(data[0][:],data[1][:],'.-b')
        plt.tight_layout()
        plt.show(block=False)     
        answer          = int(raw_input("1 to continue 0 to repeat: "))
'''
plt.clf()

################################### print ########################################
print 'QTc90' + ' = ' + str(QTc_90) +'\t'+'QTc50' + ' = ' + str(QTc_50) +'\n'
print QT_90
# sys_arg = str(sys.argv[1])

### the main results
stream = str(sys.argv[1]) + '\n' + 'QTc90,70,50,20 = ' +'\t'+ str(QTc_90)+ '\t' + str(QTc_70) + '\t' + str(QTc_50) + '\t' + str(QTc_20) + '\n'
stream += 'avgQT90,70,50,20 = ' +'\t'+ str(av_QT90)+'\t'+str(av_QT70)+'\t'+str(av_QT50)+'\t'+str(av_QT20)+'\n'
stream += 'Frequency(BPM)=' + str(BPM)+'\n'

### extras
stream1 =str(sys.argv[1]) + '\n' + 'QT90,70,50,20 = '+ '\n'
for i in range(len(QT_90)):
        stream1 += str(QT_90[i]) + '\t'
stream1 += '\n'
for i in range(len(QT_70)):
        stream1 += str(QT_70[i]) + '\t'
stream1 += '\n'
for i in range(len(QT_50)):
        stream1 += str(QT_50[i]) + '\t'
stream1 += '\n'
for i in range(len(QT_20)):
        stream1 += str(QT_20[i]) + '\t'
stream1 += '\n' + 'RR = '

for i in range(len(RR)):                 # adds RR to the stream
        stream1 += str(RR[i]) + '\t'        
stream1 = stream1 + '\n' + 'peaks = '
for i in range(len(peak_h)):                # adds RR to the stream
        stream1 = stream1 + str(peak_h[i]) + '\t'
stream1 = stream1 + '\n'
#stream = stream + str(slope) + '\n'
##################################Plotting#########################################
plot1 = plt.plot(data[0][:],data[1][:])
plot2 = plt.plot(data[0][:],baseline)
#plot10= plt.plot(data[0][:]+.5,data_f(data[0][:])+.01,'--b')
plot3 = plt.plot(data[0][peak_id],data[1][peak_id],'ro')

plot5 = plt.plot(p_x_90rl[0], data_f(p_x_90rl[0]),'go')
plot6 = plt.plot(p_x_70rl[0], data_f(p_x_70rl[0]),'bo')
plot7 = plt.plot(p_x_50rl[0], data_f(p_x_50rl[0]),'go')
plot8 = plt.plot(p_x_20rl[0], data_f(p_x_20rl[0]),'yo')
plot9 = plt.plot(p_x_50rl[1], data_f(p_x_50rl[0]),'go')
plot10 = plt.plot(p_x_90rl[1], data_f(p_x_90rl[0]),'go')
plot11 = plt.plot(p_x_70rl[1], data_f(p_x_70rl[0]),'bo')
plot12 = plt.plot(p_x_20rl[1], data_f(p_x_20rl[0]),'yo')

plt.savefig(str(sys.argv[1])+'.png')
plt.show()
with open('APD.txt', 'a') as file:
    file.write(stream)
with open('extra.txt', 'a') as file:
    file.write(stream1)
print 'Done'

