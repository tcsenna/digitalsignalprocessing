#! /bin/sh

datafile=modelinho.data

trimodel xmin=-2 zmin=0 xmax=12.0 zmax=2.0 \
         1 xedge=-2,0,2,4,6,8,10,12 \
           zedge=0,0,0,0,0,0,0,0 \
           sedge=0,0,0,0,0,0,0,0 \
         2 xedge=-2,0,2,4,6,8,10,12 \
           zedge=0.3,0.32,0.3,0.6,0.2,0.25,0.25,0.25 \
           sedge=0,0,0,0,0,0,0,0 \
         3 xedge=-2,0,2,4,6,8,10,12 \
           zedge=0.8,0.8,1.0,1.3,0.5,0.7,1.0,1.0  \
           sedge=0,0,0,0,0,0,0,0 \
         4 xedge=-2,0,2,4,6,8,10,12 \
           zedge=1.5,1.5,1.6,1.9,1.0,1.2,1.7,1.9 \
           sedge=0,0,0,0,0,0,0,0 \
         5 xedge=1.9,2.0,2.1 \
           zedge=0.4,0.36,0.4 \
           sedge=0,0,0 \
         6 xedge=1.9,2.0,2.1 \
           zedge=0.4,0.44,0.4 \
           sedge=0,0,0 \
         7 xedge=-2,0,2,4,6,8,10,12 \
           zedge=2,2,2,2,2,2,2,2 \
           sedge=0,0,0,0,0,0,0,0 \
           sfill=1,0.2,0,0,0.44,0,0 \
           sfill=1,0.5,0,0,0.15,0,0 \
           sfill=1,1.0,0,0,0.12,0,0 \
           sfill=1,1.9,0,0,0.10,0,0 \
           sfill=2,0.4,0,0,0.10,0,0 \
           kedge=1,2,3,4,5,6,7 \
           >$datafile

# Create a PS display of the model
spsplot <$datafile >modelothonia.ps \
        title="Modelinho Thonia" \
        labelz="Profundidade (km)" labelx=""\
	gedge=0.5 gtri=2.0 \
        wbox=6.0 hbox=2.0 &

exit
