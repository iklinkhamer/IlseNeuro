function current_time = Now()
time = clock;
current_time = sprintf("%d_%d_%d-%d_%d_%d", time(1), time(2), time(3), time(4), time(5), int64(time(6)));
end