pro PICOLA_mem, NPART, NGRID, BUFFER, INTSIZE, FLOATSIZE, MEM_LIMIT, MAXPROC, outfile1

; =========================================================================
; Set colors
red = reform([255, 0, 0], 1, 3)
green = reform([0, 185, 0], 1, 3)
blue = reform([0, 0, 255], 1, 3)
yellow = reform([255, 205, 0], 1, 3)
purple = reform([255, 0, 255], 1, 3)
black = reform([0, 0, 0],1 , 3)
tvlct, red, 1
tvlct, green, 2
tvlct, blue, 3
tvlct, yellow, 4
tvlct, purple, 5
tvlct, black, 6

NPART = double(NPART)
NGRID = double(NGRID)
BUFFER = double(BUFFER)
INTSIZE = double(INTSIZE)
MEM_LIMIT = double(MEM_LIMIT)
FLOATSIZE = double(FLOATSIZE)

MEM = dblarr(MAXPROC-1)
NPROC = dindgen(MAXPROC-1)+1.0;

fac_1 = 9.0*FLOATSIZE*NGRID                                                        ; 2LPT grids
fac_2 = INTSIZE*NGRID*NGRID                                                        ; 2LPT grids
fac_3 = NPART*NPART*NPART*6.0*(FLOATSIZE/2.0)                                      ; ZA and 2LPT displacements
fac_4 = NPART*NPART*NPART*BUFFER*(12.0*(FLOATSIZE/2.0)+INTSIZE)                    ; Full particle structure (including whether or not particle is inside lightcone) and ZA and 2LPT displacements
fac_5 = NPART*NPART*NPART*(BUFFER+2.0*(BUFFER-1.0))*(12.0*(FLOATSIZE/2.0)+INTSIZE) ; Particle structure including moving particles
fac_6 = 4.0*FLOATSIZE*NGRID                                                        ; Density and force grids (allocated outside the loop and so are present for moving particles and calculating displacements)
fac_7 = FLOATSIZE*NGRID*(NGRID+2.0)                                                ; Extra slice to transfer the density
fac_8 = NPART*NPART*NPART*3.0*(FLOATSIZE/2.0)                                     ; Arrays to hold the particle displacements (we allocate and deallocate these each timestep

MEM_2LPT = (fac_1/NPROC)*(NGRID*NGRID+NGRID*(2.0+NPROC)+2.0*NPROC) + fac_2
MEM_DISP = (fac_3/NPROC) + ((2.0*fac_1)/(3.0*NPROC))*(NGRID*NGRID+NGRID*(2.0+NPROC)+2.0*NPROC)
MEM_INIT = (fac_3/NPROC) + (fac_4/NPROC)
MEM_MOVE = (fac_5/NPROC)
MEM_DENS = (fac_4/NPROC) + (fac_6/NPROC)*(NGRID*NGRID+NGRID*(2.0+NPROC)+2.0*NPROC) + fac_7
MEM_NBODY = (fac_4/NPROC) + ((3.0*fac_6)/(4.0*NPROC))*(NGRID*NGRID+NGRID*(2.0+NPROC)+2.0*NPROC) + (fac_8/NPROC)

for i = 0, MAXPROC-2 do MEM[i] = max([MEM_2LPT[i],MEM_DISP[i],MEM_INIT[i],MEM_MOVE[i],MEM_DENS[i],MEM_NBODY[i]])

MEM /= (1024.0*1024.0*1024.0)
MEM_2LPT /= (1024.0*1024.0*1024.0)
MEM_DISP /= (1024.0*1024.0*1024.0)
MEM_INIT /= (1024.0*1024.0*1024.0)
MEM_MOVE /= (1024.0*1024.0*1024.0)
MEM_DENS /= (1024.0*1024.0*1024.0)
MEM_NBODY /= (1024.0*1024.0*1024.0)

index = where(MEM lt MEM_LIMIT)
outputstring = 'Minimum Number of Processors = ' + strtrim(index[0]+1, 2)

print, MEM

set_plot, 'ps'
device,filename=outfile1,/color

plot, NPROC, MEM, /nodata, background=6, xrange=[0.0, MAXPROC], yrange=[min([MEM_2LPT[MAXPROC-2],MEM_DISP[MAXPROC-2],MEM_INIT[MAXPROC-2],MEM_MOVE[MAXPROC-2],MEM_DENS[MAXPROC-2],MEM_NBODY[MAXPROC-2]])-0.2, MEM_LIMIT*2], xtitle = 'Number of Processors', ytitle='Memory Per Processor (GByte)', charthick=3, thick=3, xthick=3, ythick=3, charsize=1.5, xstyle=1, ystyle=1
oplot, NPROC, MEM_2LPT, color=1, thick=3
oplot, NPROC, MEM_DISP, color=2, thick=3
oplot, NPROC, MEM_INIT, color=3, thick=3
oplot, NPROC, MEM_MOVE, color=4, thick=3
oplot, NPROC, MEM_DENS, color=5, thick=3
oplot, NPROC, MEM_NBODY, color=6, thick=3
oplot, NPROC, dblarr(MAXPROC-1)+MEM_LIMIT, color = 1, thick = 3, linestyle = 2
xyouts, 3900, 2500, outputstring, charsize = 1.0, charthick = 3, /device

items = ['2LPT - Fields', '2LPT - Displacements', 'Initialising Particle Data', 'Moving Particles', 'Force and Density Grids', 'N-Body Displacements']
lines = [0, 0, 0, 0, 0, 0]
color = [1, 2, 3, 4, 5, 6]
al_legend, items, linestyle = lines, colors = color, number = 2, delimiter = ':', /top, /right_legend, charthick = 3, thick = 3

device, /close
set_plot, 'x'

end  


