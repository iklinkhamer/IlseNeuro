%% I.K. 7-9-23
function [path, chSel, formatSpec] = getMousePath(mname, session, whose)

switch whose
    case 'Mick'
        switch mname
            case 'One'
                chSel_type = 64;
                formatSpec_type = 1;
                switch session
                    case 1
                        path ='/home/mick/Desktop/10008-08F_2022-05-23_11-05-23_1/Record/';
                    case 2
                        disp('Second session does not exist');
                end
            case 'Two'
                chSel_type = 64;
                formatSpec_type = 1;
                switch session
                    case 1
                        path ='/home/mick/Desktop/10008-08F_2022-06-03_12-37-12_1/Record/';
                    case 2
                        disp('Second session does not exist');
                end
            case 'Three'
                chSel_type = 32;
                formatSpec_type = 2;
                switch session
                    case 1
                        path ='/home/mick/Desktop/MI13134-01_2022-12-14_12-17-28_1/Record/';
                    case 2
                        disp('Second session does not exist');
                end
        end


    case 'Ilse'
        switch mname
            case 'Apollo'
                chSel_type = 32;
                formatSpec_type = 2;
                switch session
                    case 1
                        path ='/home/mick/Desktop/spikeSortedUnits/Apollo/Session_1/';
                    case 2
                        path ='/home/mick/Desktop/spikeSortedUnits/Apollo/Session_2/';
                end
            case 'Bacchus'
                chSel_type = 32;
                formatSpec_type = 2;
                switch session
                    case 1
                        path ='/home/mick/Desktop/spikeSortedUnits/Bacchus/Session_1/';
                    case 2
                        path ='/home/mick/Desktop/spikeSortedUnits/Bacchus/Session_2/';
                end
            case 'Epona'
                switch session
                    case 1
                        chSel_type = 32;
                        formatSpec_type = 2;
                        path ='/home/mick/Desktop/spikeSortedUnits/Epona/Session_1/';
                    case 2
                        chSel_type = 64;
                        formatSpec_type = 1;
                        path ='/home/mick/Desktop/spikeSortedUnits/Epona/Session_2/';
                end
            case 'Fortuna'
                switch session
                    case 1
                        chSel_type = 64;
                        formatSpec_type = 1;
                        path ='/home/mick/Desktop/spikeSortedUnits/Fortuna/Session_1/';
                    case 2
                        chSel_type = 32;
                        formatSpec_type = 2;
                        path ='/home/mick/Desktop/spikeSortedUnits/Fortuna/Session_2/';
                end
            case 'test'
                chSel_type = 32;
                formatSpec_type = 2;
                switch session
                    case 1
                        path ='/home/mick/Desktop/spikeSortedUnits/02170-01/01/';
                    case 2
                        path ='/home/mick/Desktop/spikeSortedUnits/02170-01/02/';
                end
            case 'GoodMick'
                chSel_type = 32;
                formatSpec_type = 2;
                switch session
                    case 1
                        path = '/home/mick/Desktop/spikeSortedUnits/Mick_good_sorted/0010201/01/';
                    case 2
                        path = '/home/mick/Desktop/spikeSortedUnits/Mick_good_sorted/0010201/02/';
                end
            case 'FortunaMick'
                chSel_type = 32;
                formatSpec_type = 2;
                switch session
                    case 1
                        path = '/home/mick/Desktop/spikeSortedUnits/Fortuna/Mick_good_sorted/MI02178-04/01/';
                    case 2
                        path = '/home/mick/Desktop/spikeSortedUnits/Fortuna/Mick_good_sorted/MI02178-04/02/';
                end
        end

        if formatSpec_type == 1
            formatSpec='%s%d_CH%d.continuous'; %filename format for .continuous data
        elseif formatSpec_type ==2
            formatSpec='%s%d_%d.continuous'; %filename format for .continuous data
        end

if chSel_type == 32
    chSel = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32];
elseif chSel_type == 64
    chSel = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64];
end

end