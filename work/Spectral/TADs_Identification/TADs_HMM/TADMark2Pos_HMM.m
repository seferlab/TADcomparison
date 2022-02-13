function posn = TADMark2Pos_HMM(TAD_mark)


Labels = ones(1,length(TAD_mark));
Labels(1) = 1;
Alarm_M = 0;

for i = 2 :  length(TAD_mark)
   
   % Contiue state   1 
   Cond_cont1 =  TAD_mark(i) == 1 & TAD_mark(i-1) == 1;
   Cond_cont2 =  TAD_mark(i) == 2 & TAD_mark(i-1) == 1;
   Cond_cont3 =  TAD_mark(i) == 2 & TAD_mark(i-1) == 2;
   Cond_cont4 =  ~Alarm_M & TAD_mark(i) == 1 & TAD_mark(i-1) == 2;
   Cond_cont5 =  TAD_mark(i) == 2 & TAD_mark(i-1) == 3;
   Cond_cont = Cond_cont1 | Cond_cont2  |  Cond_cont3 | Cond_cont4 | Cond_cont5;
   
   % Alarm state    2
   Cond_alarm1 =  TAD_mark(i) == 3 & TAD_mark(i-1) == 1;
   Cond_alarm2 =  TAD_mark(i) == 3 & TAD_mark(i-1) == 2;
   Cond_alarm3 =  TAD_mark(i) == 3 & TAD_mark(i-1) == 3;
   Cond_alarm = Cond_alarm1 | Cond_alarm2 | Cond_alarm3;
   
   % Ends state & Start1
   Cond_end1 = TAD_mark(i) == 1 & TAD_mark(i-1) == 3;
   Cond_end2 = Alarm_M & (TAD_mark(i) == 1 & TAD_mark(i-1) == 2);
   Cond_end = Cond_end1 | Cond_end2;
   
   if Cond_cont
       Labels(i) = 2;
       continue;
   end
   if Cond_alarm
       Labels(i) = 2;
       Alarm_P = i;
       Alarm_M = 1;
       continue
   end
   if Cond_end
       Labels(Alarm_P) = 3;
       Labels(i) = 1;
       Alarm_M = 0;
   end
end
Labels(end) = 3;

posn = find(Labels==1);
posn(end+1) = length(Labels);
%figure,imagesc(ones(10,1)*Labels)


end