%Ian Norris Extra Credit Assignment
%AerE 161

Day = input('Please enter the day: ');
Month = input('Please enter the month: ');
Year = input('Please enter the year (YYYY): ');
fprintf('You entered: %g/%g/%g\n',Month,Day,Year)
%Initiates values for day, month, and year.

D = Day;
Y = mod(Year,100);
C = (Year - Y)/100;
%Changes format of the date to the Roman Calendar.

if Month < 3
    M = Month + 10;
    if Y == 0
        Y = 99;
        C = C - 1;
    else
        Y = Y - 1;
    end
else
    M = Month - 2;
end
%Continues changing the date to the Roman Calendar System.

A = ((13 * M) - 1)/5;
W = mod((D + fix(A) + Y + fix(Y/4) + fix(C/4) - (2*C) + 777),7);
%Calculates the value of the weekday based on the inputed date.

if W == 0
    disp('Your weekday is : Sunday')
elseif W == 1
    disp('Your weekday is : Monday')
elseif W == 2
    disp('Your weekday is : Tuesday')
elseif W == 3
    disp('Your weekday is : Wednesday')
elseif W == 4
    disp('Your weekday is : Thursday')
elseif W == 5
    disp('Your weekday is : Friday')
else
    disp('Your weekday is : Saturday')
end
%Displays the actual weekday based on the previous calculations.
