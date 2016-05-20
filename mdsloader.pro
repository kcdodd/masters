;loads the data from the files made by labview into mdsplus
; reads mdsloader_status.txt file every 5 seconds to see when it should begin loading.
; status = 0: wait another 5 seconds
; status = -1: stop waiting and close loader program
; status > 0: load data from files with shot number = status
; status file is modified by labview daq vi

pro loadfile, filename, shot
	
	print, 'Reading file: '+filename

	get_lun, lun
	openr,lun,filename

	mdsopen,'MACHINE',shot,status = stat0
	print,'mdsopen  ',shot,stat0

	line=''

	while (not EOF(lun)) do begin

		node_name = ''
		type = ''
		num_dims = ''

		readf,lun,node_name

		if (strlen(node_name) eq 0) then BREAK

		readf,lun,type
		type = fix(type,type=3)
		readf,lun,num_dims
		num_dims = fix(num_dims,type=3)

		newdata = 0

		if (type eq 0) then begin

			if (num_dims eq 1) then begin
				n = ''
				readf,lun,n
				n = fix(n,type=3)

				if (n gt 0) then begin

					data = fltarr(1, n)

					readf,lun,data

					newdata = 1
				endif
			endif

			if (num_dims eq 2) then begin
				n = ''
				readf,lun,n
				n = fix(n,type=3)	

				m = ''
				readf,lun,m
				m = fix(m,type=3)

				if ((n gt 0) and (m gt 0)) then begin

					data = fltarr(m, n)

					readf,lun,data

					newdata = 1
				endif
			endif
		endif

		if (type eq 1) then begin
				n = ''
				readf,lun,n
				n = fix(n,type=3)

				if (n gt 0) then begin

					;dataarr = StrArr(n)

					;readf,lun,dataarr
					data=''

					readf,lun,data

					for i = 0, n - 1 do begin
						;data = data+dataarr[i]+'\r\n'
					endfor

					;data = 'broken'

					newdata = 1
				endif
		end


		if newdata then begin
			print,'storing '+node_name
			mdsput,node_name,'$',data, status=st3
		endif

	endwhile


	mdsclose
	close,lun
	free_lun, lun
	
	return
end

pro mdsloader

	;initialize the status file
	get_lun,lun
	statusfile = './mdsloader_status.txt'
	openw,lun,statusfile 
	status = 0
	printf,lun,status
	close,lun
	free_lun,lun
	print,'status file created.'
	print,'Waiting for data...'


	; start waiting

	while (status eq 0) do begin

		wait, 5

		status=''
		get_lun,lun
		openr,lun,statusfile 
		readf,lun,status
		close,lun
		free_lun,lun

		status = fix(status,type=3)

		if (status gt 0) then begin

			wait, 1

			file_spec = './'+strcompress(string(status),/REMOVE_ALL)+'_spec.txt';

			loadfile,file_spec,status

			file_pres = './'+strcompress(string(status),/REMOVE_ALL)+'_pres.txt';

			loadfile,file_pres,status

			print,'data uploaded for shot:', status

			status = 0
			get_lun,lun
			openw,lun,statusfile 
			printf,lun,status
			close,lun
			free_lun,lun

			print,'Waiting for data...'

		endif

	endwhile

	print, 'daq program terminated: closing loader'

	return

end
