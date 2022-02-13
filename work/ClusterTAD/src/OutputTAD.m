Fname = strcat(Dir,algorithm,nn,'_domain.txt');
fileID = fopen(Fname,'w');

List = [];
for i = 1:length(newB)
    List =[ List , newB(i,1), newB(i,1)*Res, newB(i,2), newB(i,2)*Res ];
end
fprintf(fileID,'%6s %12s %6s %12s\n','from.id','from.cord','to.id','to.cord');

fprintf(fileID,'%6d %12d %6d %12d\n',  List);
fclose(fileID);