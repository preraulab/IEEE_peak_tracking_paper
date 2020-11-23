function my_date_str = myDateStr(min_time_unit)

clock_vec = clock;
for ii = 1:min(min_time_unit,length(clock_vec))
    if clock_vec(ii) < 10
        curr_unit = ['0' num2str(clock_vec(ii))];
    else
        curr_unit = num2str(clock_vec(ii));
    end
    if ii == 1
        my_date_str = curr_unit;
    elseif ii < 6
        my_date_str = [my_date_str curr_unit];
    else
        my_date_str = [my_date_str curr_unit(1:2)];
    end
end