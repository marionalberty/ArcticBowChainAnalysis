function Bdata = vGridBchain_plotDep(Bdata)
% Vertically grid bowchian data for 
% Restructure data to preserve original
v_0 = Bdata;
clear Bdata
Bdata.v_0 = v_0;
% Rename ungridded z & t matrix
z_ug = Bdata.v_0.z;
t_ug = Bdata.v_0.t;
% Make z grid
z1 = min([0 round(mean(z_ug(1,:)))]);
z2 = round(mean(z_ug(end,:)));
Bdata.z = (z1:-1:z2)';
% Save dnum, lat, lon to gridded version
Bdata.dn = Bdata.v_0.dn;
Bdata.lat = Bdata.v_0.lat;
Bdata.lon = Bdata.v_0.lon;
% Interp onto regular grid
Bdata.t = nan(numel(Bdata.z),numel(Bdata.dn));
for i_z = 1:numel(Bdata.dn)
  Bdata.t(:,i_z) = interp1(z_ug(:,i_z),t_ug(:,i_z),Bdata.z);
end
Bdata.z = -1*Bdata.z;