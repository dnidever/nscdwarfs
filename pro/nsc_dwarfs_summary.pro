pro nsc_dwarfs_summary,version

; Create a summary file and remove duplicates and make HTML pages

NSC_ROOTDIRS,dldir,mssdir,localdir,host
if n_elements(version) eq 0 then version = 'v3'
nside = 64
radeg = 180.0d0 / !dpi
dir = dldir+'users/dnidever/nsc/dwarfs/'+version+'/'
radius = '1.0'

npix = 49152L
pix = lindgen(npix)

;files = file_search(dir+strtrim(pix/1000,2)+'/'+strtrim(pix,2)+'_'+radius+'_'+version+'_peaks.fits')
;files = file_search(dir+strtrim(pix/1000,2)+'/'+strtrim(pix,2)+'/'+strtrim(pix,2)+'_'+radius+'_'+version+'_peaks.fits')
files = dir+strtrim(pix/1000,2)+'/'+strtrim(pix,2)+'/'+strtrim(pix,2)+'_'+radius+'_'+version+'_peaks.fits'
;test = file_test(files)
;stop
;files = file_search(dir+'*/*_peaks.fits')
;info = file_info(files)
;g = where(info.size gt 5760,ng)

;all = mrdfits(dir+'nsc_dwarf_peaks.'+version+'.ned.fits',1)
all = mrdfits(dir+'nsc_dwarf_peaks.'+version+'.labels.fits',1)
goto,mkhtml

undefine,all
for i=0,npix-1 do begin
  print,i+1,' ',files[i]
  info = file_info(files[i])
  if info.exists eq 1 then begin
    str = mrdfits(files[i],1)
    nstr = n_elements(str)
    add_tag,str,'field','',str
    add_tag,str,'mapfile','',str
    add_tag,str,'cmdfile','',str
    add_tag,str,'healpix',0L,str
    base = file_basename(files[i],'_peaks.fits')
    ;dum = strsplit(base,'_',/extract)
    ;ra1 = double(dum[0])
    ;dec1 = double(dum[1])
    ;theta = (90-dec1)/radeg
    ;phi = ra1/radeg
    theta = (90-str.dec)/radeg
    phi = str.ra/radeg
    ang2pix_ring,nside,theta,phi,pix
    str.field = base
    str.mapfile = base+'_skymap_peaks.png'
    str.cmdfile = base+'_cmd'+strtrim(lindgen(nstr)+1,2)+'.png'
    str.healpix = pix[0]
    if n_elements(all) eq 0 then begin
      schema = str
      struct_assign,{dum:''},schema
      all = replicate(schema,npix)
    endif
    all[i] = str
    ;push,all,str
  endif
endfor

stop

; only keep ones that succeeded
gd = where(all.id ne '',ngd)
all = all[gd]
; 26530 of 

; Only keeps peaks that are inside their healpix boundary
gd = where(all.hpix eq all.healpix,ngd)
; 14641

; Merge duplicate detections of a single object due to sub-clumps?

; Add galactic coordinates
add_tag,all,'glon',0.0d0,all
add_tag,all,'glat',0.0d0,all
glactc,all.ra,all.dec,2000.0,glon,glat,1,/deg
all.glon = glon
all.glat = glat

;mwrfits,all,dir+'nsc_dwarf_peaks.'+version+'.fits',/create




stop

mkhtml:

; make html page
hdir = dir+'html/'
nall = n_elements(all)
nperpage = 1000
npages = nall/nperpage+1
cnt1 = 0LL
for p=0,npages-1 do begin
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
  push,lines,'<tr><th>NUM</th><th>RA</th><th>DEC</th><th>PEAK</th><th>LABEL</th><th>MAP</th><th>CMD</th></tr>'
  cnt = n_elements(lines)

  push,lines,strarr(9*nperpage)
  for i=0,nperpage-1 do begin
    ;push,lines,strarr(8*n_elements(all))
    ;for i=0,n_elements(all)-1 do begin
    if cnt1 lt nall then begin
      if i mod 10 eq 0 then print,i
      lines[cnt] = '<tr>'
      lines[cnt+1] = '<td>'+strtrim(cnt1+1,2)+'</td>'
      lines[cnt+2] = '<td>'+stringize(all[cnt1].ra,ndec=5)+'</td>'
      lines[cnt+3] = '<td>'+stringize(all[cnt1].dec,ndec=5)+'</td>'
      lines[cnt+4] = '<td>'+stringize(all[cnt1].peak_value,ndec=3)+'</td>'
      lines[cnt+5] = '<td>'+strtrim(all[cnt1].label_type,2)+'</td>'
      lines[cnt+6] = '<td><img src="'+all[cnt1].mapfile+'"></td>'
      lines[cnt+7] = '<td><img src="'+all[cnt1].cmdfile+'"></td>'
      lines[cnt+8] = '</tr>'
      cnt += 9
    endif
    cnt1++
  endfor
  push,lines,'</body>'
  push,lines,'</html>'
  writeline,hdir+'nsc_dwarf_peaks'+strtrim(p+1,2)+'.html',lines
  ;writeline,'nsc_dwarf_peaks.html',lines

endfor

; Make index page
undefine,lines
push,lines,'<html>'
push,lines,'<head>'
push,lines,'<title>'
push,lines,'</title>'
push,lines,'</head>'
push,lines,'<body>'
push,lines,'<h3>NSC Dwarf Search<h3>'
for i=0,npages-1 do push,lines,'<a href="nsc_dwarf_peaks'+strtrim(i+1,2)+'.html">Page '+strtrim(i+1,2)+'</a><br>'
push,lines,'</body>'
push,lines,'</html>'
writeline,hdir+'index.html',lines

stop

file_copy,'nsc_dwarf_peaks1.html','html1'
ui=uniq(all[0:99].mapfile,sort(all[0:99].mapfile))
file_copy,all[ui].mapfile,hdir,/over
file_copy,all[0:99].cmdfile,hdir,/over

;file_copy,'nsc_dwarf_peaks.html','html'
ui=uniq(all.mapfile,sort(all.mapfile))
file_copy,dir+strtrim(all[ui].hpix/1000,2)+'/'+strtrim(all[ui].hpix,2)+'/'+all[ui].mapfile,hdir,/over
;file_copy,dir+all[ui].mapfile,hdir,/over
ui=uniq(all.cmdfile,sort(all.cmdfile))
file_copy,dir+strtrim(all[ui].hpix/1000,2)+'/'+strtrim(all[ui].hpix,2)+'/'+all[ui].cmdfile,hdir,/over
;file_copy,all.cmdfile,hdir,/over


stop

end
