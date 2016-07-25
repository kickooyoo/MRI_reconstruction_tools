function dcf = get_dcf(Nro, Nspokes)
	dcfRow=ones(Nro,1);
	for ii=1:Nro
		dcfRow(ii)=abs(Nro/2-(ii-0.5));
	end
	dcfScl=norm(dcfRow);
	dcfRow=1/dcfScl*dcfRow;
	dcf=repmat(dcfRow,1,Nspokes);
end	
