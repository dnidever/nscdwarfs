pro nsc_dwarfs_summary

; Create a summary file and remove duplicates and make HTML pages

NSC_ROOTDIRS,dldir,mssdir,localdir,host
version = 'v2'
nside = 64
radeg = 180.0d0 / !dpi
dir = dldir+'users/dnidever/nsc/dwarfs/'+version+'/'
radius = '1.0'

pix = lindgen(49152)

files = file_search(dir+strtrim(pix/1000,2)+'/'+strtrim(pix,2)+'_'+radius+'_'+version+'_peaks.fits')
test = file_test(files)
stop
;files = file_search(dir+'*/*_peaks.fits')
info = file_info(files)
g = where(info.size gt 5760,ng)

undefine,all
for i=0,ng-1 do begin
  print,i+1,' ',files[g[i]]
  str = mrdfits(files[g[i]],1)
  nstr = n_elements(str)
  add_tag,str,'field','',str
  add_tag,str,'mapfile','',str
  add_tag,str,'cmdfile','',str
  add_tag,str,'healpix',0L,str
  base = file_basename(files[g[i]],'_peaks.fits')
  dum = strsplit(base,'_',/extract)
  ra1 = double(dum[0])
  dec1 = double(dum[1])
  theta = (90-dec1)/radeg
  phi = ra1/radeg
  ang2pix_ring,nside,theta,phi,pix
  str.field = base
  str.mapfile = base+'_skymap_peaks.png'
  str.cmdfile = base+'_cmd'+strtrim(lindgen(nstr)+1,2)+'.png'
  str.healpix = pix[0]
  push,all,str
endfor

;mwrfits,all,'nsc_dwarf_peaks.fits',/create
;9019 peaks.

; Only keeps peaks that are inside their healpix boundary

stop

; make html page
undefine,lines
push,lines,'<html>'
push,lines,'<head>'
push,lines,'<title>'
push,lines,'</title>'
push,lines,'</head>'
push,lines,'<body>'
push,lines,'<h3>NSC Dwarf Search<h3>'
push,lines,'<hr>'
push,lines,'<table border=1>'
push,lines,'<tr><th>NUM</th><th>RA</th><th>DEC</th><th>PEAK</th><th>MAP</th><th>CMD</th></tr>'
cnt = n_elements(lines)
n = 100
push,lines,strarr(8*n)
for i=0,n-1 do begin
;push,lines,strarr(8*n_elements(all))
;for i=0,n_elements(all)-1 do begin
if i mod 10 eq 0 then print,i
lines[cnt] = '<tr>'
lines[cnt+1] = '<td>'+strtrim(i+1,2)+'</td>'
lines[cnt+2] = '<td>'+stringize(all[i].ra,ndec=5)+'</td>'
lines[cnt+3] = '<td>'+stringize(all[i].dec,ndec=5)+'</td>'
lines[cnt+4] = '<td>'+stringize(all[i].peak_value,ndec=3)+'</td>'
lines[cnt+5] = '<td><img src="'+all[i].mapfile+'"></td>'
lines[cnt+6] = '<td><img src="'+all[i].cmdfile+'"></td>'
lines[cnt+7] = '</tr>'
cnt += 8
endfor
push,lines,'</body>'
push,lines,'</html>'
writeline,'nsc_dwarf_peaks1.html',lines
;writeline,'nsc_dwarf_peaks.html',lines

file_copy,'nsc_dwarf_peaks1.html','html1'
ui=uniq(all[0:99].mapfile,sort(all[0:99].mapfile))
file_copy,all[ui].mapfile,'html1',/over
file_copy,all[0:99].cmdfile,'html1',/over

file_copy,'nsc_dwarf_peaks.html','html'
ui=uniq(all.mapfile,sort(all.mapfile))
file_copy,all[ui].mapfile,'html',/over
file_copy,all.cmdfile,'html',/over


stop

end
