'''
Plot functions for the Brake Squeal Project

'''

import matplotlib.pyplot as plt
def plot_eigs_cover(obj,la):

        x1 = obj.target[0]
	x2 = obj.target[1]
	y1 = obj.target[2]
	y2 = obj.target[3]
	
	# calculate center of the obj.target rectangular region
        tau = complex((x1+x2)/2,(y1+y2)/2)
        
        radius=max(abs(la-tau))
         
	fig = plt.gcf()
	circle=plt.Circle((tau.real,tau.imag),radius,color='b',fill=False)
	react = plt.Rectangle((x1,y1),abs(x2-x1),abs(y2-y1),color='black',fill=False)
	ax = plt.gca()
	ax.cla() # clear things for fresh plot
	ax.set_xlim((tau.real-1.1*radius,tau.real+1.1*radius))
	ax.set_ylim((tau.imag-1.1*radius,tau.imag+1.1*radius))
	plt.axhline(0, color='blue')
	plt.axvline(0, color='blue')
	x = la.real
	y = la.imag
	for i in range(0,len(x)):
		if x[i] >= 0:
			ax.plot(x[i],y[i],'o',color='red')
		else:
			ax.plot(x[i],y[i],'o',color='green')	
	
	ax.plot((tau.real),(tau.imag),'+',color='r')
	fig.gca().add_artist(circle)
	fig.gca().add_artist(react)
	fig.savefig(obj.output_path+'plotcover.png')
	return radius;

def plot_eigs_transition(obj,la):
        
        x1 = obj.target[0]
	x2 = obj.target[1]
	y1 = obj.target[2]
	y2 = obj.target[3]
	
	# calculate center of the obj.target rectangular region
        tau = complex((x1+x2)/2,(y1+y2)/2)
        
        radius=max(abs(la-tau))
	fig = plt.gcf()
	circle=plt.Circle((tau.real,tau.imag),radius,color='b',fill=False)
	react = plt.Rectangle((x1,y1),abs(x2-x1),abs(y2-y1),color='black',fill=False)
	ax = plt.gca()
	ax.cla() # clear things for fresh plot
	ax.set_xlim((-200,200))
	ax.set_ylim((y1-100,y2+100))
	plt.axhline(0, color='blue')
	plt.axvline(0, color='blue')
	x = la.real
	y = la.imag
	for i in range(0,len(x)):
		if x[i] >= 0:
			ax.plot(x[i],y[i],'.',color='red')
		else:
			ax.plot(x[i],y[i],'.',color='green')	
	
	ax.plot((tau.real),(tau.imag),'+',color='r')
	#fig.gca().add_artist(circle)
	#fig.gca().add_artist(react)
	fig.savefig(obj.output_path+'plottransition.png')
	return radius;
