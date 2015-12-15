%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Riveted Joint Graphing Code         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = dlmread('Specimen_RawData.txt');
time = zeros(3776,1);
extension = zeros(3776,1);
load = zeros(3776,1);
for y = 1:3776
    for z = 1:3
        if(z==1)
            time(y) = x(y,z);  %Seconds
        elseif(z==2)
            extension(y) = x(y,z); %Millimeters
        else
            load(y) = x(y,z);   %Newtons
        end
    end
end
            
length = 203.2;   %Millimeters
Area = .508*2*6.35; %Millimeters^2
stress = load./(Area/(1000^2));  %Newtons/Meters^2
strain = extension./length;   %Unitless  (Millimeters/Millimeters)

figure(1)
plot(stress,strain)
title('Stress vs. Strain')
xlabel('Stress (N/m^2)')
ylabel('Strain (mm/mm)')

figure(2)
plot(time,stress)
title('Time vs. Stress')
xlabel('Time (seconds)')
ylabel('Stress (N/m^2)')

figure(3)
plot(time,strain)
title('Time vs. Strain')
xlabel('Time (seconds)')
ylabel('Strain (mm/mm)')