
import matplotlib.pyplot as plt
import numpy as np
import math  

####################### functions ###############
def centerPlot(cellCenterX, cellCenterY) : 

    plt.scatter(cellCenterX, cellCenterY, color="k")

    plt.xlabel('X (micro meter)')
    plt.ylabel('Y (micro meter)')
    plt.title('locations of cell cetner')
    plt.grid(True)
    plt.savefig("test.png")
    plt.show() 


####################### input parameter data  ###############
nPouchCells=65# 29  
nPeripCells=11 #6
nBcCells=5 # each side

lumen=2 #2 

wPouch=2.5 #2 
hPouchA=25 #10   height of the paouch cells at the front of cell center if we rotate CCW
hPouchB=25 #10   height of the paouch cells at the front of cell center if we rotate CCW

bNodePouch=9 #18 #14  #number of basal node
aNodePouch=9 #18 #14 #number of apical nodes
lNodePouchA=aNodePouch*10    #5
lNodePouchB=aNodePouch*10    #5

wPerip=16.4 #10 +0.92783  
hPerip=4  #2 

bNodePerip=4*14 #4*28  #number of basal node
aNodePerip=4*14 #4*28 #number of apical nodes
lNodePerip=14 #28   #14


domainMinX=-1.394 
cellsGapX=0.329  
domainMinY=14.3  
domainMaxY=domainMinY+hPouchA+lumen+hPerip  

bNodeBc=14 #30 #12  #number of basal node  ## needs to be even number
aNodeBc=14 #30 #12 #number of apical nodes ## needs to be even number
lNodeBc=14 #30 #12 ## needs to be even number

 
 # 100 ## for one side nPouchCells=29  
correctionFactorECMBCY=-0.25

buffer_ECM=0.525 
rBc=2.25 ## radius of boundary cells initialized with circular shape
dhECM=0.1
enlargeR=rBc+buffer_ECM


################## calculation of coordinate of cell centers ######################
cellCenterPouchX=[domainMinX+cellsGapX+wPouch/2] 
cellCenterPouchY=[domainMinY+hPouchA/2] 
cellCenterPouchTheta=[0]
for x in range(nPouchCells-1):
    cellCenterPouchX.append ( cellCenterPouchX[x]+wPouch + cellsGapX)
    cellCenterPouchY.append ( cellCenterPouchY[x]) 
    cellCenterPouchTheta.append (0)

cellCenterPeripX =[domainMinX+cellsGapX+wPerip/2]  
cellCenterPeripY=[domainMinY+hPouchA+lumen+hPerip/2] 
cellCenterPeripTheta=[np.pi]
for x in range (nPeripCells-1) :
    cellCenterPeripX.append ( cellCenterPeripX[x]+wPerip + cellsGapX)   
    cellCenterPeripY.append ( cellCenterPeripY[x] )  
    cellCenterPeripTheta.append (np.pi)

radiusBc=(cellCenterPeripY[0]-cellCenterPouchY[0])/2
centerRightBcX=cellCenterPouchX[nPouchCells-1] + wPouch/2 + cellsGapX
centerRightBcY=0.5* ( cellCenterPeripY[0]+cellCenterPouchY[0] ) 


cellCenterBcX=[0 for x in range (nBcCells)]
cellCenterBcY=[0 for x in range (nBcCells)]
cellCenterBcTheta=[0 for x in range (nBcCells)]

correctFactorA=1.2
correctFactorB=np.pi*0.1

for x in range (nBcCells) : 
    cellCenterBcX[x]= centerRightBcX+ radiusBc*np.sin ( (x+1)*np.pi*correctFactorA/(nBcCells+1) -correctFactorB)
    cellCenterBcY[x]= centerRightBcY- radiusBc*np.cos ( (x+1)*np.pi*correctFactorA/(nBcCells+1) -correctFactorB)
    cellCenterBcTheta[x]=(x+1)*np.pi*correctFactorA/(nBcCells+1)


centerLeftBcX=cellCenterPouchX[0] - wPouch/2-cellsGapX
centerLeftBcY=centerRightBcY


for x in range (nBcCells) : 
    cellCenterBcX.append ( centerLeftBcX- radiusBc*np.sin ((x+1)*np.pi*correctFactorA/(nBcCells+1)-correctFactorB) )
    cellCenterBcY.append ( centerLeftBcY+ radiusBc*np.cos ((x+1)*np.pi*correctFactorA/(nBcCells+1)-correctFactorB) )
    cellCenterBcTheta.append ( (x+1)*np.pi*correctFactorA/(nBcCells+1) +np.pi )

######################### unify and sort cell centers locations ############
cellCenterX=[]
cellCenterY=[]
cellType=[]
cellCenterTheta=[]

for x in range (nPouchCells) : 
    cellCenterX.append (cellCenterPouchX[x] )  
    cellCenterY.append (cellCenterPouchY[x] )
    cellType.append ('pouch')  
    cellCenterTheta.append (0)

for x in range (nBcCells) : 
    cellCenterX.append ( cellCenterBcX[x] )  
    cellCenterY.append ( cellCenterBcY[x] )
    cellType.append ('bc')
    cellCenterTheta.append (cellCenterBcTheta[x])

for x in range (1,nPeripCells+1) :
    cellCenterX.append ( cellCenterPeripX[nPeripCells-x] )  
    cellCenterY.append ( cellCenterPeripY[nPeripCells-x] ) 
    cellType.append ('peri') 
    cellCenterTheta.append (np.pi)

for x in range (nBcCells,2*nBcCells) : 
    cellCenterX.append ( cellCenterBcX[x] )  
    cellCenterY.append ( cellCenterBcY[x] )
    cellType.append ('bc')
    cellCenterTheta.append (cellCenterBcTheta[x])

####################### plot cell centers locations #####################
centerPlot(cellCenterX, cellCenterY) 


####################### Write cell centers locations as an output file #####################

fileM = open('coordinate_Cell16.txt','w') 

numCells= nPouchCells + nPeripCells + 2*nBcCells 
fileM.write('{}\n'.format(numCells))
for i in range ( numCells ) :
    #fileM.write(str(i) +"," +str (xC[i][j])   )  #,yC[i][j],typeC[i][j]) 
    fileM.write('{0:.3f}'.format(cellCenterX[i]))
    fileM.write(' ')
    fileM.write('{0:.3f}'.format(cellCenterY[i]))
    fileM.write(' ')

    fileM.write('{0:.3f}'.format(0.000))
    fileM.write(' ')

    fileM.write('{}\n'.format(cellType[i]))


fileM.close()



########################### start finding the membrane nodes ##################################################


xC=[]  ## for k =0 to create the list
yC=[]  ## for k =0 to create the list
typeC=[]
dpp=[]


############ finding membrane nodes location of pouch cells #######################
for k in range ( nPouchCells ) :
    lastpoint=0 
    if k==0 :
        lNodePouchB=18 #36 #14
        lNodePouchA=54 #108 #42
        aNodePouch=20  #40
        bNodePouch=20  #40
        hPouchB=5  #2
        hPouchA=15 #6
            
    elif  k==(nPouchCells-1) :
        lNodePouchB=54 #108 #42
        lNodePouchA=18 #36 #14
        aNodePouch=20 #40
        bNodePouch=20 #40
        hPouchB=15   #15 #6
        hPouchA=5  #2 
    elif k==1 :
        lNodePouchB=54 #108 #42
        lNodePouchA=90 #180
        aNodePouch=20 #40
        bNodePouch=20 #40
        hPouchB=15 #6
        hPouchA=25 #10  
    elif k==(nPouchCells-2) :
        lNodePouchB=90 #180
        lNodePouchA=54 #108 # 42
        aNodePouch=20 #40
        bNodePouch=20 #40
        hPouchB=25 #10
        hPouchA=15 #6
    else:
        lNodePouchB=90 #180
        lNodePouchA=90 #180 # 42
        aNodePouch=9 #18
        bNodePouch=9 #18
        hPouchB=25 #10
        hPouchA=25 #6

       
    dh=0.5*abs (hPouchA-hPouchB)
    cellTotalNodePouch=lNodePouchB+lNodePouchA+bNodePouch+aNodePouch
    xC.append ( [0 for x in range ( cellTotalNodePouch)] )
    yC.append ( [0 for x in range ( cellTotalNodePouch)] )
    dpp.append ( [0 for x in range ( cellTotalNodePouch)] )
    typeC.append ( ['notAssigned1' for x in range ( cellTotalNodePouch)] ) 
    
    # from center to top equal to H/2
    for i  in range (int(lNodePouchA/2))  :
        j=lastpoint+i 
        typeC[k][j]='lateralA'
        xC[k][j]=cellCenterX[k]+wPouch/2 
        yC[k][j]=cellCenterY[k]+i/(lNodePouchA/2)*hPouchA/2  

    # on apical side
    lastpoint=lastpoint+int (lNodePouchA/2)  
    for i  in range (aNodePouch) :
        j=lastpoint+i 
        typeC[k][j]='apical1'
        xC[k][j]=cellCenterX[k]+wPouch/2- i/aNodePouch*wPouch 
        if k==0 or k==1 :
            yC[k][j]=cellCenterY[k]+hPouchA/2 - dh*(i)/aNodePouch 
        elif k==(nPouchCells-2) or k==(nPouchCells-1) :
            yC[k][j]=cellCenterY[k]+hPouchA/2 + dh*i/aNodePouch  
        else:
            yC[k][j]=cellCenterY[k]+hPouchA/2 


    lastpoint=lastpoint+ aNodePouch  

    for i  in range (lNodePouchB) :
        j=lastpoint+i 
        typeC[k][j]='lateralB'
        xC[k][j]=cellCenterX[k]-wPouch/2.0 
        yC[k][j]=cellCenterY[k]+hPouchB/2.0  - i/lNodePouchB*hPouchB 


    lastpoint=lastpoint+ lNodePouchB 
    for i  in range (bNodePouch) :
        j=lastpoint+i 
        typeC[k][j]='basal1'
        xC[k][j]=cellCenterX[k]-wPouch/2+ i/bNodePouch*wPouch

        if k==0 or k==1 :
            yC[k][j]=cellCenterY[k]-hPouchB/2 -dh*i/bNodePouch
        elif k==(nPouchCells-2) or k==(nPouchCells-1) :
            yC[k][j]=cellCenterY[k]-hPouchB/2 +dh*(i)/bNodePouch
        else:
            yC[k][j]=cellCenterY[k]-hPouchB/2 


    lastpoint=lastpoint+ bNodePouch 
    for i  in range (int(lNodePouchA/2))  :
        j=lastpoint+i 
        typeC[k][j]='lateralA'
        xC[k][j]=cellCenterX[k]+wPouch/2 
        yC[k][j]=cellCenterY[k]-hPouchA/2+ i/(lNodePouchA/2)*hPouchA/2  

############ finding membrane nodes location of right hand side BC cells #######################

for k in range ( nPouchCells,nPouchCells+int(nBcCells) ) :
    lastpoint=0

    cellTotalNodeBc=2*lNodeBc+bNodeBc+aNodeBc
    typeC.append ( ['notAssigned1' for x in range ( cellTotalNodeBc)] ) 
    xC.append ( [0 for x in range ( cellTotalNodeBc)] )
    yC.append ( [0 for x in range ( cellTotalNodeBc)] ) 
    dpp.append ( [0 for x in range ( cellTotalNodeBc)] )

    # from center to top equal to H/2
    for i  in range (int(lNodeBc/2))  :
        j=lastpoint+i 
        typeC[k][j]='lateralA'
        angle=(i+1)/(lNodeBc/2)*np.pi*0.25  +cellCenterTheta[k]
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(angle)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(angle)


    lastpoint=lastpoint+int (lNodeBc/2) 
    lastAngle=angle 
    for i  in range (int(aNodeBc))  :
        j=lastpoint+i 
        typeC[k][j]='apical1'
        angle=(i+1)/aNodeBc*np.pi*0.5+lastAngle
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(angle)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(angle)

    lastpoint=lastpoint+int (aNodeBc)  
    lastAngle=angle
    for i  in range (int(lNodeBc))  :
        j=lastpoint+i
        typeC[k][j]='lateralB'
        angle=(i+1)/lNodeBc*np.pi*0.5+lastAngle
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(angle)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(angle)


    lastpoint=lastpoint+int (lNodeBc) 
    lastAngle=angle 
    for i  in range (int(bNodeBc))  :
        j=lastpoint+i 
        typeC[k][j]='basal1'
        angle=(i+1)/bNodeBc*np.pi*0.5+lastAngle
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(angle)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(angle)


    lastpoint=lastpoint+int (bNodeBc)
    lastAngle=angle   
    for i  in range (int(lNodeBc/2))  :
        j=lastpoint+i 
        typeC[k][j]='lateralA'
        angle=(i+1)/(lNodeBc/2)*np.pi*0.25+lastAngle
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(angle)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(angle)


############ finding membrane nodes location of peripodial cells #######################

for k in range ( nPouchCells+int(nBcCells),nPouchCells+int(nBcCells)+nPeripCells ) :

    lastpoint=0
    cellTotalNodePerip=2*lNodePerip+bNodePerip+aNodePerip
    typeC.append ( ['notAssigned1' for x in range ( cellTotalNodePerip)] )
    xC.append ( [0 for x in range ( cellTotalNodePerip)] )
    yC.append ( [0 for x in range ( cellTotalNodePerip)] ) 
    dpp.append ( [0 for x in range ( cellTotalNodePerip)] )
    # from center to top equal to H/2
    for i  in range (int(lNodePerip/2))  :
        j=lastpoint+i
        typeC[k][j]='lateralA'
        xC[k][j]=cellCenterX[k]-wPerip/2 
        yC[k][j]=cellCenterY[k]-i/(lNodePerip/2)*hPerip/2  


    lastpoint=lastpoint+int (lNodePerip/2)  
    for i  in range (aNodePerip) :
        j=lastpoint+i 
        typeC[k][j]='apical1'
        xC[k][j]=cellCenterX[k]-wPerip/2+ i/aNodePerip*wPerip
        yC[k][j]=cellCenterY[k]-hPerip/2  


    lastpoint=lastpoint+ aNodePerip  
    for i  in range (lNodePerip) :
        j=lastpoint+i
        typeC[k][j]='lateralB' 
        xC[k][j]=cellCenterX[k]+wPerip/2 
        yC[k][j]=cellCenterY[k]-hPerip/2  + i/lNodePerip*hPerip 
    

    lastpoint=lastpoint+ lNodePerip 
    for i  in range (bNodePerip) :
        j=lastpoint+i 
        typeC[k][j]='basal1' 
        xC[k][j]=cellCenterX[k]+wPerip/2- i/bNodePerip*wPerip
        yC[k][j]=cellCenterY[k]+hPerip/2  


    lastpoint=lastpoint+ bNodePerip 
    for i  in range (int(lNodePerip/2))  :
        j=lastpoint+i
        typeC[k][j]='lateralA' 
        xC[k][j]=cellCenterX[k]-wPerip/2 
        yC[k][j]=cellCenterY[k]+hPerip/2- i/(lNodePerip/2)*hPerip/2  

############ finding location of membrane nodes of left hand side BC cells #######################

for k in range ( nPouchCells+nPeripCells+int(nBcCells),nPouchCells+nPeripCells+int(2*nBcCells) ) :

    lastpoint=0
    cellTotalNodeBc=2*lNodeBc+bNodeBc+aNodeBc
    typeC.append ( ['notAssigned1' for x in range ( cellTotalNodeBc)] )
    xC.append ( [0 for x in range ( cellTotalNodeBc)] )
    yC.append ( [0 for x in range ( cellTotalNodeBc)] ) 
    dpp.append ( [0 for x in range ( cellTotalNodeBc)] )
    # from center to top equal to H/2
    for i  in range (int(lNodeBc/2))  :
        j=lastpoint+i 
        typeC[k][j]='lateralA' 
        angle=(i+1)/(lNodeBc/2)*np.pi*0.25 +cellCenterTheta[k]
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(angle)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(angle)


    lastpoint=lastpoint+int (lNodeBc/2)  
    lastAngle=angle
    for i  in range (int(aNodeBc))  :
        j=lastpoint+i 
        typeC[k][j]='apical1'
        angle=(i+1)/aNodeBc*np.pi*0.5+lastAngle 
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(angle)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(angle)


    lastpoint=lastpoint+int (aNodeBc)  
    lastAngle=angle
    for i  in range (int(lNodeBc))  :
        j=lastpoint+i 
        typeC[k][j]='lateralB' 
        angle=(i+1)/lNodeBc*np.pi*0.5+lastAngle
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(angle)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(angle)


    lastpoint=lastpoint+int (lNodeBc)  
    lastAngle=angle
    for i  in range (int(bNodeBc))  :
        j=lastpoint+i 
        typeC[k][j]='basal1'
        angle=(i+1)/bNodeBc*np.pi*0.5+lastAngle
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(angle)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(angle)


    lastpoint=lastpoint+int (bNodeBc)  
    lastAngle=angle
    for i  in range (int(lNodeBc/2))  :
        j=lastpoint+i 
        typeC[k][j]='lateralA' 
        angle=(i+1)/(lNodeBc/2)*np.pi*0.25+lastAngle
        xC[k][j]=cellCenterX[k]+ rBc*np.cos(angle)
        yC[k][j]=cellCenterY[k]+ rBc*np.sin(angle)




## calculate DPP level for each cells ##################################
numNode=[0 for x in range (numCells)]
#num_Node=[ 0 for x in range (numCells) ]
totalNode=0 
for k in range ( numCells ) :
    numNode [k]= len (xC[k]) 
    totalNode=totalNode+numNode [k]


counter=0
sumCoordX=0
for i in range ( numCells ) :
    for j in range (numNode [i]): # number of membrane nodes for that cell
        sumCoordX=sumCoordX+ xC[i][j]
        counter=counter+1 

tissueCenterX=sumCoordX/counter

maxEachCell=[0 for x in range (numCells)]
for i in range ( numCells ) : 
        maxEachCell[i]=np.amax (xC[i])

maxXC=np.amax ( maxEachCell)
        


for i in range ( numCells ) :
    for j in range (numNode [i]): # number of membrane nodes for that cell
       dpp[i][j]=max ( 1.0 - abs (tissueCenterX-xC[i][j])/abs (tissueCenterX-maxXC),0.0 )

######################## plot membrane nodes locations #########################
for k in range ( numCells ) :   
    
    plt.scatter(xC[k], yC[k], c=dpp[k], vmin=0, vmax=1)
    plt.gca().set_aspect('equal', adjustable='box')
#plt.xlabel('X (micro meter)')
#plt.ylabel('Y (micro meter)')
#plt.title('locations of cell cetner')
#plt.grid(True)
#plt.savefig("test.png")
#plt.show() 

#centerPlot(xC[0], yC[0]) 
######################## write membrane nodes location and type #########################



fileM = open('coordinate_Membrane7.txt','w') 
 
#fileM.write('cellID,x coordinate, y coordinate, nodeType\n') 
for i in range ( numCells ) :
    for j in range (numNode [i]): # number of membrane nodes for that cell
        #fileM.write(str(i) +"," +str (xC[i][j])   )  #,yC[i][j],typeC[i][j]) 
        fileM.write('{}'.format(i))
        fileM.write(' ')
        fileM.write('{0:.4f}'.format(xC[i][j]))
        fileM.write(' ')
        fileM.write('{0:.4f}'.format(yC[i][j]))
        fileM.write(' ')
        fileM.write('{0:.4f}'.format(dpp[i][j]))
        fileM.write(' ')
        fileM.write('{}\n'.format(typeC[i][j]))
       
        

fileM.close() 


########################### start finding the ECM nodes location and type ##################################################
xECM=[]
yECM=[]
eCMType=[]
numECMnodesCount=0



## ###########ECM nodes location below the non-rectangular pouch cells at the Left ################



maxY=0 # small number
for i in range (len(yC[0])) :   
    if yC[0][i]>maxY and typeC[0][i]=='basal1' :
        maxY=yC[0][i]

x1=min (xC[0])   
y1=maxY  ## of basal nodes in cell with id zero 
x2=(max(xC[1]))
y2=(min (yC[1]))

m=(y2-y1)/(x2-x1)
b=y1-m*x1


xECM.append(x1-0.15) 
yECM.append(y1-buffer_ECM) 
eCMType.append ('excm')
i=0
while xECM[numECMnodesCount]<max (xC[1]):
    i=i+1
    xTmp=x1+i*dhECM/math.sqrt(1+m*m)
    xECM.append(xTmp-0.15)      #/max (xC[62])-min (xC[2])) 
    yECM.append(m*xTmp+b-buffer_ECM)
    eCMType.append ('excm')
    numECMnodesCount+=1


################## pouch cells  uniform in size


i=-1
while xECM[numECMnodesCount]<max (xC[62]):
    i=i+1
    xECM.append(min (xC[2]) +i*dhECM)        #/max (xC[62])-min (xC[2])) 
    yECM.append(min (yC[2])-buffer_ECM)
    eCMType.append ('excm')
    numECMnodesCount+=1


# ###########ECM nodes location below the non-rectangular pouch cells at the Right ################



maxY=0  # small number
for i in range (len(yC[64])) :   
    if yC[64][i]>maxY and typeC[64][i]=='basal1' :
        maxY=yC[0][i]

x1=min (xC[63])   
y1=min (yC[63])
x2=(max(xC[64]))
y2=maxY # for the basal nodes

m=(y2-y1)/(x2-x1)
b=y1-m*x1

i=-1
while xECM[numECMnodesCount]<max (xC[64]):
    i=i+1
    xTmp=x1+i*dhECM/math.sqrt(1+m*m)
    xECM.append(xTmp+0.1)    
    yECM.append(m*xTmp+b-buffer_ECM)
    eCMType.append ('excm')
    numECMnodesCount+=1

################# calculation of coordinates of ECM nodes which are neighbor with right hand side BC cells ####################
numECMBCNodes=int ( np.pi *(radiusBc+enlargeR)/dhECM)
for k in range ( numECMBCNodes) :
    xECM.append(centerRightBcX +(radiusBc+enlargeR)*np.sin ( (k+1)*np.pi/(numECMBCNodes+1) ))
    yECM.append(centerRightBcY+correctionFactorECMBCY- (radiusBc+enlargeR)*np.cos ( (k+1)*np.pi/(numECMBCNodes+1) ))
    eCMType.append ('bc2')
    numECMnodesCount+=1

################# calculation of coordinates of peripodial side ECM ####################


i=-1 # local counter
while xECM[numECMnodesCount]>min (xC[80]):
    i=i+1
    xECM.append(max (xC[70]) -i*dhECM)       
    yECM.append(max (yC[70])+buffer_ECM)
    eCMType.append ('perip')
    numECMnodesCount+=1

################# calculation of coordinates of ECM nodes which are neighbor with left hand side BC cells ####################
for k in range ( numECMBCNodes) :
    xECM.append(  centerLeftBcX- (radiusBc+enlargeR)*np.sin ((k+1)*np.pi/(numECMBCNodes+1)) )
    yECM.append(  centerLeftBcY+correctionFactorECMBCY+ (radiusBc+enlargeR)*np.cos ((k+1)*np.pi/(numECMBCNodes+1)) ) 
    eCMType.append ('bc2')
    numECMnodesCount+=1



########################### plotting ECM nodes ##################################################
plt.scatter(xECM, yECM, color="k")
plt.show()                


########################### writing as an output file ECM nodes locations and type ##################################################
fileM = open('coordinate_ECM16.txt','w') 
 
#fileM.write('cellID,x coordinate, y coordinate, nodeType\n')
fileM.write('{}\n'.format(numECMnodesCount)) 
for k in range ( numECMnodesCount ) :
        #fileM.write(str(i) +"," +str (xC[i][j])   )  #,yC[i][j],typeC[i][j]) 
        fileM.write('{0:.4f}'.format(xECM[k]))
        fileM.write(' ')
        fileM.write('{0:.4f}'.format(yECM[k]))
        fileM.write(' ')
        fileM.write('{}\n'.format(eCMType[k]))


fileM.close() 











