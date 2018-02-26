#Does many runs of tV_Diagonalize for dfferent M and N
#NOTE: Currently written for half-filling case only
import subprocess
import numpy as np

program = 'tV_main.jl'
N = 2          #Max number of particles
xA = '2'              #Partition size
xAout = 'l2'
case = 'F'
siteMax = '1'
code = 'EOP'
bc = '--pbc'
bcOut = 'to be determined'
uMin = '-1.522878745'  #log_10(0.030000000019365026)
uMax = '2'   #log_10(100)
#uMin = '1'  #log_10(0.030000000019365026)
#uMax = '1'   #log_10(100)
uNum = '80' #Number of points to be calculated
alpha = '2'        #Renyi Index
negLabel = 'yes'   #Type 'yes' if doing negative U run. Otherwise, 'no'.

#Paste the right part of the following at the end of subprocess.call() if desired:
#'--probs'    #Generate file with probabilities respect to subsystem size
#'--neg'      #Do run using negatives of the specified values

#Add NEG label to negative U runs
if negLabel=='yes':
    extension = 'NEG.dat'
    
else:
    extension = '.dat'

for n in range(N,N+1,9873):

    #Change between antiperiodic and periodic BC label
    if n % 2 == 0:
        bcOut = 'A'
    else:
        bcOut =  'P'

    m = n*2
    m = str(m)
    n = str(n)

    output = code+bcOut+m+case+n+xAout+'a'+alpha+extension
    
    #Comment/Uncomment desired option:

    #1) Log Scale Call
    subprocess.call(['julia', program, '--u-log', bc, '--u-num', uNum, '--u-min', uMin, '--u-max', uMax,'--site-max', siteMax, m, n, '--ee', xA, '--out', output,'--alpha', alpha,'--neg'])

    #2) Linear Scale Call
    #subprocess.call(['julia', program, bc, '--u-num', uNum, '--u-min', uMin, '--u-max', uMax,'--site-max', siteMax, m, n, '--ee', xA, '--out', output,'--alpha', alpha])

