%% plot percentage of tuned units per 50ms bins for motor imagery (SMG, PMV and S1 - figure 1C/D) and speech (SMG  - Figure 4D)

%% Important: run code while being in folder 'grasp_and_speech_decoding'
addpath(genpath(pwd)); %add folder to search path 

%%
clc
clear all 
close all

%%
% parameters
phase_names = {'ITI', 'Cue', 'Delay', 'Action'};
number_phases = length(phase_names); 

figures_to_plot = {'1C', '1D', '4D'};
Fr = 20; %frequency rate of firing rate

for n_f = 1:length(figures_to_plot)
    
%load data corresponding to figure that is being plotted

    if contains(figures_to_plot{n_f}, '1')
        
        %load Go or NoGo motor imagery data
        if strcmp(figures_to_plot{n_f}, '1C')
           motor_imagery_data = load([pwd '\Data\UnitTuning50msBins\GoData_MotorImagery']);
           plot_title = 'Go trials';

        elseif strcmp(figures_to_plot{n_f}, '1D')
            motor_imagery_data = load([pwd '\Data\UnitTuning50msBins\NoGoData_MotorImagery']);     
            plot_title = 'No-Go trials';

        end 
            Data1 = motor_imagery_data.data1;
            Data2 = motor_imagery_data.data2;
            Data3 = motor_imagery_data.data3;
            total_number_units = motor_imagery_data.num_channels_total;
            labels ={'SMG', 'PMV', 'S1X'};
    end 
        
    if strcmp(figures_to_plot{n_f}, '4D')
        %load speech data
        speech_data = load([pwd '\Data\UnitTuning50msBins\Speech']);
        %load motor imagery data
        motor_imagery_data = load([pwd '\Data\UnitTuning50msBins\GoData_MotorImagery']);
        
        Data1 = motor_imagery_data.data1; %all motor imagery data
        Data2 = speech_data.data1; % grasp speech data
        Data3 = speech_data.data2; % color speech data
        total_number_units = [motor_imagery_data.num_channels_total(1),speech_data.num_channels_total(1), speech_data.num_channels_total(2)];

        plot_title = 'SMG';
        labels = {'Motor Imagery', 'Spoken Grasps', 'Spoken Colors'};

    end 
    
    time_phase_labels = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2;...
        2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;...
        2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;...
        4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;...
        4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4];
    time_phase_labels_adapted = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;0;0;0;0;0;...
                                2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;...
                                0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;...
                                0;0;0;0;0;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;0;0;0;0;0;0;0;0;0;0;0;0;...
                                0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
    
    time = (1:length(Data1))/Fr;
    phase_labels(1)=1/Fr;
    phase_labels(2:4) = (find(diff(time_phase_labels))+1)/Fr;
    
    %figure variables
    xtitle_font_size =13;
    title_font_size = 14;
    font_size = 15;
    colors_tasks = util.get_color_rgb_codes(labels);
    y_lim =70;
    %generate labels including number of units
    labels2 = arrayfun(@(x) [labels{x} ' (N = ' num2str(total_number_units(x)) ')'], 1:length(total_number_units), 'UniformOutput', false);

    %plot each task phase in a different subplot
    fig = figure('units','normalized','outerposition',[0 0 0.45 0.45]);
    for k = 1:length(phase_names)
        if k  == 1
            subplot(1,6,1)
        elseif k == 2
            subplot(1,6,2:3)
        elseif k == 3
            subplot(1,6,4)
        elseif k == 4
            subplot(1,6,5:6)
        end 

        sub_time_idx = find(time_phase_labels == k);
        sub_time = time(sub_time_idx);
        plot(sub_time, Data1(sub_time_idx),'Color',colors_tasks{1}, 'LineWidth', 1.5); 
        hold on 
        plot(sub_time, Data2(sub_time_idx),'Color',colors_tasks{2},'LineWidth', 1.5); 
        hold on 
        plot(sub_time, Data3(sub_time_idx), 'Color',colors_tasks{3}, 'LineWidth', 1.5);
        ylim([0, y_lim+2])
        title(phase_names{k}, 'FontSize', font_size)
        xlim([sub_time(1), sub_time(end) ])

        hold on
        if k == 2|| k == 4
            plot(sub_time(find(time_phase_labels_adapted(sub_time_idx))), ones(size(sub_time(find(time_phase_labels_adapted(sub_time_idx)))))*70, 'LineWidth', 1, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15,'Color', [100 100 100]/255 )
        end 

        %plot legend only for 4th panel
        if k == 4
            legend(labels2, 'FontSize', font_size)
        end 

        %change the xticks for the smaller plots 
        try        
            if k == 1 || k == 3
                xt = xticks;
                xticks(xt([2,4]))
            end
        catch
            %do nothing
        end 

        %remove y ticks for all the plots that are not the first subplot
        if k ~=1
            set(gca,'YTick', [])
        end 
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',12) 

        %xlabel, ylabel and title of figure
        han=axes(fig,'visible','off'); 
        han.Title.Visible='on';
        han.XLabel.Visible='on';
        han.YLabel.Visible='on';
        ylabel(han,'Tuned units [%]', 'FontSize', font_size);
        xlabel(han,'Time [s]', 'FontSize', font_size);
        t = title(han,plot_title,'FontSize', font_size);
    end 

end 