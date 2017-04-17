%% robust k-space trajectory mapping
%  supports 2D or 3D
%

%% clean slate
clear all; close all; clc;

%% get path to this mfile
mfile_fullpath = mfilename('fullpath');
[folder_parent,~,~] = fileparts(mfile_fullpath);

%% add paths
addpath( sprintf('%s/thirdparty/sdc3_nrz_11aug/', folder_parent) );
addpath( sprintf('%s/thirdparty/grid3_dct_11aug/', folder_parent) );

%% set folder paths
data_input_folder = sprintf('%s/data_input', folder_parent);
data_output_folder = sprintf('%s/data_output', folder_parent);
png_folder = sprintf('%s/png', folder_parent);

%% trajectory type list
trajectory_type_list = {'spiral','radial','multivane','cartesian','epi'};

%% common settings
field_of_view_mm = 256;
acquired_voxel_size_mm = 1;
slice_distance_to_isocenter_mm = 5;
        
for trajectory_idx = 1:numel(trajectory_type_list),

    %% set trajectory type
    trajectory_type = trajectory_type_list{trajectory_idx};
    
    %% set data_input_folder and filename according to trajectory_type
    switch trajectory_type,
        case 'cartesian'
            filename = '20160913_125634_PB_TFE_cartesian_PN_XYZ_0_0_0_offset_5mm.mat';
            acquired_voxel_size_mm = 1;

        case 'spiral'
            filename = '20160913_130002_PB_TFE_spiral_PN_XYZ_0_0_0_offset_5mm.mat';
            acquired_voxel_size_mm = 1;

        case 'radial'
            filename = '20160913_124914_PB_TFE_radial_PN_XYZ_0_0_0_offset_5mm.mat';
            acquired_voxel_size_mm = 1;

        case 'multivane'
            filename = '20160913_125301_PB_TFE_multivane_PN_XYZ_0_0_0_offset_5mm.mat';
            acquired_voxel_size_mm = 1;

        case 'epi'
            filename = '20160913_130044_PB_EPI_PN_XYZ_0_0_0_offset_5mm_thick_2mm.mat';
            acquired_voxel_size_mm = 4;
    end

    %% load .MAT file (created from LAB/RAW/SIN using ReadPhilipstoMAT.py
    loaded = load( sprintf('%s/%s', data_input_folder, filename) );

    % converted .MAT data is nReps x nKy x nKx
    [nReps nKy nKx] = size(loaded.data);
    ksp_all = complex(zeros([nKx nKy nReps]));
    for rep = 1:nReps,
        for ky = 1:nKy,
            ksp_all(:,ky,rep) = loaded.data(rep,ky,:);
        end
    end
    nKz = 1;

    %% Calculate number of dynamic averages
    % x,y,z +1,-1 has multiple of 13 reps 
    nAvgs = nReps / 13;

    %% average k-space
    ksp = complex(zeros([nKx nKy 12]));
    ksp_encoded = complex(zeros([nKx nKy 1]));
    for avg = 1:nAvgs,    
        ksp = ksp + ksp_all(:,:, [1:12] + 13*(avg-1) )/nAvgs;
        ksp_encoded = ksp_encoded + ksp_all(:,:, 13 + 13*(avg-1) )/nAvgs;
    end

    %% display original encoded k-space
    hfig = figure(1); clf;
    set(hfig,'Position',[600 900 1000 1000]);
    subplot(2,2,1); imagesc(log(abs(ksp_encoded))); axis image; colormap(gray); title('encoded k-space (abs)');
    subplot(2,2,2); imagesc(angle(ksp_encoded)); axis image; colormap(gray); title('encoded k-space (angle)');
    subplot(2,2,3); imagesc(real(ksp_encoded)); axis image; colormap(gray); title('encoded k-space (real)');
    subplot(2,2,4); imagesc(imag(ksp_encoded)); axis image; colormap(gray); title('encoded k-space (imag)');

    %% save figure to PNG file
    [pathstr, prefix, ext] = fileparts(filename);
    imwrite(frame2im(getframe(hfig)), sprintf('%s/%s_ksp_encoded.png', png_folder, prefix));

    %% trajecory specific k-space massaging
    switch trajectory_type

        case {'cartesian'}
            % phase chopping of encoded k-space
            for ky=2:2:nKy,
                ksp_encoded(:,ky) = ksp_encoded(:,ky) * exp(-1i*pi);
            end

        case {'epi'}
            % flip frequency encode direction
            ksp_encoded = flipud(ksp_encoded);
            for rep=1:12,
                ksp(:,:,rep) = flipud(ksp(:,:,rep));
            end

    end

    switch trajectory_type    
        case {'multivane','propeller'}
            bladewidth = 16;
        case {'cartesian','epi'}
            bladewidth = nKy;
    end

    switch trajectory_type,
        case {'radial','spiral'}

            ref_v = zeros(nKy,1);
            for ky = 1:nKy,
               [maxval, maxidx] = max(abs(ksp_encoded(:,ky)));
               ref_v(ky) = maxidx;
            end
        
            %% plot ksp_abs
            hfig = figure(2); clf;
            [pathstr, prefix, ext] = fileparts(filename);
            set(hfig,'Position',[600 900 1000 1000]);
            plot(abs(ksp(:,:))); 
            axis([1 size(ksp,1) 0 max(abs(ksp(:)))]);
            grid on; axis square;
            ht = title( sprintf('Duyn readout magnitude (%s)', prefix) );
            set(ht,'interpreter','none');
            xlabel('sample no.'); ylabel('signal magnitude (a.u.)');
            imwrite(frame2im(getframe(hfig)), sprintf('%s/%s_Duyn_readout_magnitude.png', png_folder, prefix));
        
            %% just keep phase with unity magnitude for ksp
            ksp = exp(-1i*angle(ksp));
        
            %% calculate k-space coordinates
            % x,y,z +1,-1 has multiple of 13 reps 

            % x
            x_diff_slice_pos = zeros(nKx,nKy);
            x_diff_slice_neg = zeros(nKx,nKy);
            for ky = 1:nKy,
                x_diff_slice_pos(:,ky) = unwrap_from_ref( angle( ksp(:,ky,1)./ksp(:,ky,2) ), ref_v(ky) );
                x_diff_slice_neg(:,ky) = unwrap_from_ref( angle( ksp(:,ky,3)./ksp(:,ky,4) ), ref_v(ky) );
            end

            kx_coords = x_diff_slice_pos - x_diff_slice_neg;

            % y
            y_diff_slice_pos = zeros(nKx,nKy);
            y_diff_slice_neg = zeros(nKx,nKy);
            for ky = 1:nKy,
                y_diff_slice_pos(:,ky) = unwrap_from_ref( angle( ksp(:,ky,5)./ksp(:,ky,6) ), ref_v(ky) );
                y_diff_slice_neg(:,ky) = unwrap_from_ref( angle( ksp(:,ky,7)./ksp(:,ky,8) ), ref_v(ky) );
            end

            ky_coords = y_diff_slice_pos - y_diff_slice_neg;       

            % z
            z_diff_slice_pos = zeros(nKx,nKy);
            z_diff_slice_neg = zeros(nKx,nKy);
            for ky = 1:nKy,
                z_diff_slice_pos(:,ky) = unwrap_from_ref( angle( ksp(:,ky,9)./ksp(:,ky,10) ), ref_v(ky) );
                z_diff_slice_neg(:,ky) = unwrap_from_ref( angle( ksp(:,ky,11)./ksp(:,ky,12) ), ref_v(ky) );
            end

            kz_coords = z_diff_slice_pos - z_diff_slice_neg;        

        case {'multivane','propeller','cartesian','epi'}
            bladecount = nKy/bladewidth;

            ksp_reformat = complex(zeros([nKx*bladewidth bladecount 12]));

            for rep = 1:12,
                for ky = 1:nKy,

                    blade_line = mod(ky-1,bladewidth) + 1;
                    blade_number = ceil(ky/bladewidth);                

                    switch trajectory_type,
                        case {'multivane','propeller'}
                            data_ksp_v = ksp(:, (blade_line-1)*bladewidth + blade_number, rep);
                        case {'cartesian','epi'}
                            data_ksp_v = ksp(:, ky, rep);
                    end

                    % flip every other line to keep endpoints adjacent in k-space
                    if mod(ky,2)==0,
                        data_ksp_v = flip(data_ksp_v);
                    end

                    ksp_reformat([1:nKx] + (blade_line - 1) * nKx, blade_number, rep) = data_ksp_v;

                end

            end

            ksp_encoded_reformat = complex(zeros([nKx*bladewidth bladecount]));
            ksp_encoded_reformat_tmp = complex(zeros([nKx bladewidth*bladecount]));
            for ky = 1:nKy,

                blade_line = mod(ky-1,bladewidth) + 1;
                blade_number = ceil(ky/bladewidth);

                switch trajectory_type,
                    case {'multivane','propeller'}
                        data_ksp_encoded_v = ksp_encoded(:, (blade_line-1)*bladewidth + blade_number);
                    case {'cartesian','epi'}
                        data_ksp_encoded_v = ksp_encoded(:, ky);
                end

                ksp_encoded_reformat_tmp([1:nKx], blade_line + (blade_number-1)*bladewidth ) = data_ksp_encoded_v;

                % flip every other line to keep endpoints adjacent in k-space
                if mod(ky,2)==0,
                    data_ksp_encoded_v = flip(data_ksp_encoded_v);
                end

                ksp_encoded_reformat([1:nKx] + (blade_line - 1) * nKx, blade_number) = data_ksp_encoded_v;

            end

            ref_v = zeros(bladecount,1);
            for blade = 1:bladecount,
               [maxval, maxidx] = max(abs(ksp_encoded_reformat(:,blade)));
               ref_v(blade) = maxidx;
            end

            %% plot ksp_abs
            hfig = figure(2); clf;
            [pathstr, prefix, ext] = fileparts(filename);
            set(hfig,'Position',[600 900 1000 1000]);
            ksp_reformat_abs_2D = abs(ksp_reformat(:,:));
            switch trajectory_type,
                case {'cartesian','multivane'}
                    plot(flipud(ksp_reformat_abs_2D));
                case {'epi'}
                    plot(ksp_reformat_abs_2D);
            end
            axis([1 size(ksp_reformat,1) 0 max(abs(ksp_reformat(:)))]);
            grid on; axis square;
            ht = title( sprintf('Duyn readout magnitude (%s)', prefix) );
            set(ht,'interpreter','none');
            xlabel('sample no.'); ylabel('signal magnitude (a.u.)');
            switch trajectory_type,
                case {'epi'}
                    legend('1','2','3','4','5','6','7','8','9','10','11','12');
            end
            imwrite(frame2im(getframe(hfig)), sprintf('%s/%s_Duyn_readout_magnitude.png', png_folder, prefix));

            %% just keep phase with unity magnitude for ksp
            ksp_reformat = exp(-1i*angle(ksp_reformat));

            %% calculate k-space coordinates
            % x,y,z +1,-1 has multiple of 13 reps

            % x
            x_diff_slice_pos = zeros(nKx*bladewidth,bladecount);
            x_diff_slice_neg = zeros(nKx*bladewidth,bladecount);
            for blade = 1:bladecount,
                x_diff_slice_pos(:,blade) = unwrap_from_ref( angle( ksp_reformat(:,blade,1)./ksp_reformat(:,blade,2) ), ref_v(blade) );
                x_diff_slice_neg(:,blade) = unwrap_from_ref( angle( ksp_reformat(:,blade,3)./ksp_reformat(:,blade,4) ), ref_v(blade) );
            end

            kx_coords = x_diff_slice_pos - x_diff_slice_neg;

            % y
            y_diff_slice_pos = zeros(nKx*bladewidth,bladecount);
            y_diff_slice_neg = zeros(nKx*bladewidth,bladecount);
            for blade = 1:bladecount,
                y_diff_slice_pos(:,blade) = unwrap_from_ref( angle( ksp_reformat(:,blade,5)./ksp_reformat(:,blade,6) ), ref_v(blade) );
                y_diff_slice_neg(:,blade) = unwrap_from_ref( angle( ksp_reformat(:,blade,7)./ksp_reformat(:,blade,8) ), ref_v(blade) );
            end

            ky_coords = y_diff_slice_pos - y_diff_slice_neg;

            % z
            z_diff_slice_pos = zeros(nKx*bladewidth,bladecount);
            z_diff_slice_neg = zeros(nKx*bladewidth,bladecount);
            for blade = 1:bladecount,
                z_diff_slice_pos(:,blade) = unwrap_from_ref( angle( ksp_reformat(:,blade,9)./ksp_reformat(:,blade,10) ), ref_v(blade) );
                z_diff_slice_neg(:,blade) = unwrap_from_ref( angle( ksp_reformat(:,blade,11)./ksp_reformat(:,blade,12) ), ref_v(blade) );
            end

            kz_coords = z_diff_slice_pos - z_diff_slice_neg;
    end

    %% calculate phase from B0 eddy currents
    phase_B0_eddy_x = (x_diff_slice_pos + x_diff_slice_neg)/4;
    phase_B0_eddy_y = (y_diff_slice_pos + y_diff_slice_neg)/4;
    phase_B0_eddy_z = (z_diff_slice_pos + z_diff_slice_neg)/4;

    %% plot phase from B0 eddy currents
    hfig = figure(22); clf;
    [pathstr, prefix, ext] = fileparts(filename);
    set(hfig,'Position',[200 200 1700 800]);
    subplot(1,3,1); plot(phase_B0_eddy_x * 180 / pi); grid on; title('phase from B0 eddy currents (X axis)');
    axis([1 size(phase_B0_eddy_x,1) -60  60]);
    xlabel('sample no.'); ylabel('phase (degrees)');
    subplot(1,3,2); plot(phase_B0_eddy_y * 180 / pi); grid on; title('phase from B0 eddy currents (Y axis)');
    axis([1 size(phase_B0_eddy_y,1) -60  60]);
    xlabel('sample no.'); ylabel('phase (degrees)');
    subplot(1,3,3); plot(phase_B0_eddy_z * 180 / pi); grid on; title('phase from B0 eddy currents (Z axis)');
    axis([1 size(phase_B0_eddy_y,1) -60  60]);
    xlabel('sample no.'); ylabel('phase (degrees)');
    imwrite(frame2im(getframe(hfig)), sprintf('%s/%s_phase_from_B0_eddy_currents.png', png_folder, prefix));

    %% save detected phase from B0 eddy currents
    [pathstr, prefix, ext] = fileparts(filename);
    save( sprintf('%s/%s_phase_B0_eddy_x_y_z.mat', data_output_folder, prefix), 'phase_B0_eddy_x', 'phase_B0_eddy_y', 'phase_B0_eddy_z');

    %% convert coordinates from phase to reciprocal distance units
    kmax_inverse_mm = (1/acquired_voxel_size_mm) / 2.0;
    kx_coords_inverse_mm = kx_coords / (2 * pi * slice_distance_to_isocenter_mm) / 4; % factor of 4 due to doubling the phase difference twice, subtraction of subtraction
    ky_coords_inverse_mm = ky_coords / (2 * pi * slice_distance_to_isocenter_mm) / 4;
    kz_coords_inverse_mm = kz_coords / (2 * pi * slice_distance_to_isocenter_mm) / 4;

    %% plot figure
    hfig = figure(3); clf;
    set(hfig,'Position',[200 200 1200 1000]);
    subplot(2,2,1);
    plot(kx_coords);
    title('kx'); xlabel('sample no.'); ylabel('k');
    axis square; grid on;

    subplot(2,2,2);
    plot(ky_coords);
    title('ky'); xlabel('sample no.'); ylabel('k');
    axis square; grid on;

    subplot(2,2,3);
    plot(kz_coords);
    title('kz'); xlabel('sample no.'); ylabel('k');
    axis square; grid on;

    subplot(2,2,4);
    plot3(kx_coords_inverse_mm, ky_coords_inverse_mm, kz_coords_inverse_mm);
    title('kx, ky'); xlabel('kx [mm^{-1}]'); ylabel('ky [mm^{-1}]'); zlabel('kz [mm^{-1}]');
    axis square; grid on;
    axis(kmax_inverse_mm*([-1 +1 -1 +1 -1 +1]*1.05));

    %% save figure to PNG file
    [pathstr, prefix, ext] = fileparts(filename);
    imwrite(frame2im(getframe(hfig)), sprintf('%s/%s_trajectory.png', png_folder, prefix));

    %% apply rotation matrix from .csv file (taken from location_matrices entry of .sin file)
    tmp = load( sprintf('%s/%s.csv', data_input_folder, prefix) );    
    R = zeros(3,3);
    R(1,:) = tmp(1,:);
    R(2,:) = tmp(2,:);
    R(3,:) = tmp(3,:);
    rotated_kx_coords_inverse_mm = zeros(size(kx_coords_inverse_mm));
    rotated_ky_coords_inverse_mm = zeros(size(ky_coords_inverse_mm));
    rotated_kz_coords_inverse_mm = zeros(size(kz_coords_inverse_mm));
    for idx = 1:numel(kx_coords_inverse_mm),
       v = R * [kx_coords_inverse_mm(idx) ; ky_coords_inverse_mm(idx) ; kz_coords_inverse_mm(idx)];
       rotated_kx_coords_inverse_mm(idx) = v(1);
       rotated_ky_coords_inverse_mm(idx) = v(2);
       rotated_kz_coords_inverse_mm(idx) = v(3);
    end

    %% plot figure
    hfig = figure(33); clf;
    set(hfig,'Position',[200 200 1200 1000]);
    subplot(2,2,1);
    plot(rotated_kx_coords_inverse_mm);
    title('rotated kx'); xlabel('sample no.'); ylabel('k mm^{-1}');
    axis square; grid on;

    subplot(2,2,2);
    plot(rotated_ky_coords_inverse_mm);
    title('rotated ky'); xlabel('sample no.'); ylabel('k mm^{-1}');
    axis square; grid on;

    subplot(2,2,3);
    plot(rotated_kz_coords_inverse_mm);
    title('kz'); xlabel('sample no.'); ylabel('k mm^{-1}');
    axis square; grid on;

    subplot(2,2,4);
    plot(rotated_kx_coords_inverse_mm, rotated_ky_coords_inverse_mm);
    title('rotated kx, rotated ky'); xlabel('kx [mm^{-1}]'); ylabel('ky [mm^{-1}]');
    axis square; grid on;
    axis(kmax_inverse_mm*([-1 +1 -1 +1]*1.05));
    hold on;
    plot([-kmax_inverse_mm -kmax_inverse_mm], [-kmax_inverse_mm +kmax_inverse_mm],'r-');
    plot([-kmax_inverse_mm +kmax_inverse_mm], [+kmax_inverse_mm +kmax_inverse_mm],'r-');
    plot([+kmax_inverse_mm +kmax_inverse_mm], [+kmax_inverse_mm -kmax_inverse_mm],'r-');
    plot([+kmax_inverse_mm -kmax_inverse_mm], [-kmax_inverse_mm -kmax_inverse_mm],'r-');

    %% save figure to PNG file
    [pathstr, prefix, ext] = fileparts(filename);
    imwrite(frame2im(getframe(hfig)), sprintf('%s/%s_trajectory_rotated.png', png_folder, prefix));

    %% reconstruct encoded k-space

    %% normalize by calculated kmax_inverse_mm
    % some values max be outside [-0.5,+0.5]
    kx_coords_normalized = rotated_kx_coords_inverse_mm / kmax_inverse_mm * 0.5;
    ky_coords_normalized = rotated_ky_coords_inverse_mm / kmax_inverse_mm * 0.5;
    %kz_coords_normalized = rotated_kz_coords_inverse_mm / kmax_inverse_mm * 0.5;
    kz_coords_normalized = 0.0 * rotated_kz_coords_inverse_mm / kmax_inverse_mm * 0.5; % force zero rotated z coords
    disp( sprintf('range (kx,ky,kz) = [%.3f, %.3f]', max([kx_coords_normalized(:) ; ky_coords_normalized(:) ; kz_coords_normalized(:)]), min([kx_coords_normalized(:) ; ky_coords_normalized(:) ; kz_coords_normalized(:)]) ) );

    %% save normalized coordinates
    [pathstr, prefix, ext] = fileparts(filename);
    save( sprintf('%s/%s_kx_ky_kz_coords_normalized.mat', data_output_folder, prefix), 'kx_coords_normalized', 'ky_coords_normalized', 'kz_coords_normalized');

    %% kx, ky coordinates
    switch trajectory_type,
        case {'radial','spiral'}
            crds = zeros(3, nKx, nKy);
        case {'multivane','propeller','cartesian','epi'}
            crds = zeros(3, nKx*bladewidth, nKy/bladewidth);
            ksp_encoded = ksp_encoded_reformat;
    end
    crds(1,:,:) = kx_coords_normalized;
    crds(2,:,:) = ky_coords_normalized;
    crds(3,:,:) = kz_coords_normalized;

    %% effective Matrix size
    effMtx  = 256;

    %% roll off kernel
    delta = [1.0, 0.0];
    k_not = [0.0, 0.0, 0.0];
    numThread = 1 % only have 1 data point
    rokern = grid3_MAT(delta',k_not',[1.0],effMtx,numThread);
    % change to complex, shift, then fft
    rokern = squeeze(rokern(1,:,:,:) + 1j*rokern(2,:,:,:));
    rokern = sum(rokern,3);
    rokern = fft2(fftshift(rokern));
    rokern = fftshift(rokern);
    rokern = abs(rokern);

    %% display rokern
    hfig = figure(4); clf;
    set(hfig,'Position',[600 900 1000 1000]);
    imagesc(rokern); axis image; colormap(gray);
    title('roll off kernel');

    %% save figure to PNG file
    [pathstr, prefix, ext] = fileparts(filename);
    imwrite(frame2im(getframe(hfig)), sprintf('%s/%s_rokern.png', png_folder, prefix));

    %% calculate density compensation
    numIter = 25;
    osf     = 2.1;
    verbose = 1;
    DCF = sdc3_MAT(crds,numIter,effMtx,verbose,osf);

    %% grid data
    numThread = 1; % or 1 for no pthread lib exec
    data = zeros(2, nKx, nKy);
    data(1,:) = real(ksp_encoded(:));
    data(2,:) = imag(ksp_encoded(:));
    gdata = grid3_MAT(data,crds,DCF,effMtx,numThread);

    %% change to complex, fft, then shift
    gdata_complex = squeeze(gdata(1,:,:,:) + 1j*gdata(2,:,:,:));

    gdata_complex = fftshift(gdata_complex,1);
    gdata_complex = fftshift(gdata_complex,2);
    gdata_complex = fftshift(gdata_complex,3);

    gdata_3D = ifftn(gdata_complex);

    gdata_3D = fftshift(gdata_3D,1);
    gdata_3D = fftshift(gdata_3D,2);
    gdata_3D = fftshift(gdata_3D,3);

    gdata_img = gdata_3D(:,:,129);
    %gdata_img = sum(gdata_3D,3);

    %% apply rolloff
    gdata_img(rokern > 0) = gdata_img(rokern > 0) ./ rokern(rokern > 0);

    %% crop according to trajectory_type
    switch trajectory_type,
        case {'cartesian'}
            nKx_crop = floor(nKx / 2 * 0.90);
            crop_idx = [1:nKx_crop] + size(gdata_img,1)/2 - nKx_crop/2;
            gdata_img = gdata_img(crop_idx,crop_idx);     
        case {'epi'}
            nKx_crop = size(loaded.data,2) + 1;
            crop_idx = [1:nKx_crop] + size(gdata_img,1)/2 - nKx_crop/2;
            gdata_img = gdata_img(crop_idx,crop_idx);
    end

    %% display k-space of gdata_img
    hfig = figure(5); clf
    set(hfig,'Position',[600 900 1000 1000]);
    imagesc( log( abs(fftshift(fft2(gdata_img))) ) ); axis image; colormap(gray); title('log( abs(fftshift(fft2(gdata_img))) )');

    %% save figure to PNG file
    [pathstr, prefix, ext] = fileparts(filename);
    imwrite(frame2im(getframe(hfig)), sprintf('%s/%s_gridded_kspace.png', png_folder, prefix));

    %% disaply gdata_img
    hfig = figure(6); clf;
    set(hfig,'Position',[600 900 1400 1100]);
    subplot(2,2,1); imagesc(abs(gdata_img)); axis image; colormap(gray); title('magnitude');
    subplot(2,2,2); imagesc(angle(gdata_img)); axis image; colormap(gray); title('phase');
    subplot(2,2,3); imagesc(real(gdata_img)); axis image; colormap(gray); title('real');
    subplot(2,2,4); imagesc(imag(gdata_img)); axis image; colormap(gray); title('imag');

    %% save figure to PNG file
    [pathstr, prefix, ext] = fileparts(filename);
    imwrite(frame2im(getframe(hfig)), sprintf('%s/%s_image.png', png_folder, prefix));

    %% reconstruct again compendating for phase from B0 eddy currents
    ksp_encoded_B0EC_comp = ksp_encoded .* exp(+1i * (phase_B0_eddy_x + phase_B0_eddy_y + phase_B0_eddy_z) ); % seems to be help

    data = zeros(2, nKx, nKy);
    data(1,:) = real(ksp_encoded_B0EC_comp(:));
    data(2,:) = imag(ksp_encoded_B0EC_comp(:));
    gdata = grid3_MAT(data,crds,DCF,effMtx,numThread);

    %% change to complex, fft, then shift
    gdata_complex = squeeze(gdata(1,:,:,:) + 1j*gdata(2,:,:,:));

    gdata_complex = fftshift(gdata_complex,1);
    gdata_complex = fftshift(gdata_complex,2);
    gdata_complex = fftshift(gdata_complex,3);

    gdata_3D = ifftn(gdata_complex);

    gdata_3D = fftshift(gdata_3D,1);
    gdata_3D = fftshift(gdata_3D,2);
    gdata_3D = fftshift(gdata_3D,3);

    gdata_img_B0EC_comp = gdata_3D(:,:,129);

    %% apply rolloff
    gdata_img_B0EC_comp(rokern > 0) = gdata_img_B0EC_comp(rokern > 0) ./ rokern(rokern > 0);

    %% crop according to trajectory_type
    switch trajectory_type,
        case {'cartesian'}
            nKx_crop = floor(nKx / 2 * 0.90);;
            crop_idx = [1:nKx_crop] + size(gdata_img_B0EC_comp,1)/2 - nKx_crop/2;
            gdata_img_B0EC_comp = gdata_img_B0EC_comp(crop_idx,crop_idx);    
        case {'epi'}
            nKx_crop = size(loaded.data,2) + 1;
            crop_idx = [1:nKx_crop] + size(gdata_img_B0EC_comp,1)/2 - nKx_crop/2;
            gdata_img_B0EC_comp = gdata_img_B0EC_comp(crop_idx,crop_idx);
    end

    mask_average = ( abs(gdata_img) + abs(gdata_img_B0EC_comp) )/2;
    mask_average = mask_average - min(mask_average(:));
    mask_average = mask_average / max(mask_average(:));
    mask = mask_average > graythresh(mask_average);
    background_mask_idx = find(~imdilate(mask,strel('disk',5),'same'));

    %% display gdata_img vs. gdata_img_B0EC_comp
    hfig = figure(7); clf;
    set(hfig,'Position',[50 50 1500 800]);
    subplot(1,3,1); imagesc(abs(gdata_img)); axis image; colormap(gray); title('No B0 Eddy Current Compensation');
    subplot(1,3,2); imagesc(abs(gdata_img_B0EC_comp)); axis image; colormap(gray); title( sprintf('With B0 Eddy Current Compensation\nBackground Energy Ratio (comp)/(not comp) = %.3f', ...
        sum(abs(gdata_img_B0EC_comp(background_mask_idx)))/sum(abs(gdata_img(background_mask_idx))) ) );
    subplot(1,3,3); imagesc((abs(gdata_img)-abs(gdata_img_B0EC_comp))); axis image; colormap(gray); title('Difference');

    %% save figure to PNG file
    [pathstr, prefix, ext] = fileparts(filename);
    imwrite(frame2im(getframe(hfig)), sprintf('%s/%s_image_with_without_B0EC_comp.png', png_folder, prefix));

    %% save PNG for ISMRM poster
    fs = 32;
    fw = 'bold';
    hfig = figure(1000); clf;
    set(hfig,'Position',[200 200 1000 1000]);
    subplot(1,1,1);
    plot(rotated_kx_coords_inverse_mm, rotated_ky_coords_inverse_mm);
    title('rotated kx, rotated ky'); xlabel('kx [mm^{-1}]'); ylabel('ky [mm^{-1}]');
    set(gca,'FontSize', fs, 'FontWeight', 'bold');
    axis square; grid on;
    axis(kmax_inverse_mm*([-1 +1 -1 +1]*1.05));
    hold on;
    plot([-kmax_inverse_mm -kmax_inverse_mm], [-kmax_inverse_mm +kmax_inverse_mm],'r-');
    plot([-kmax_inverse_mm +kmax_inverse_mm], [+kmax_inverse_mm +kmax_inverse_mm],'r-');
    plot([+kmax_inverse_mm +kmax_inverse_mm], [+kmax_inverse_mm -kmax_inverse_mm],'r-');
    plot([+kmax_inverse_mm -kmax_inverse_mm], [-kmax_inverse_mm -kmax_inverse_mm],'r-');
    imwrite(frame2im(getframe(hfig)), sprintf('%s/%s_rotated_trajectory_ISMRM_poster.png', png_folder, prefix));

    %% image with B0 eddy current compensation
    img = abs(gdata_img_B0EC_comp);
    img = img - min(img(:));
    img = img / max(img(:));
    imwrite(img, sprintf('%s/%s_image_with_B0EC_comp_ISMRM_poster.png', png_folder, prefix));
    
end % end trajectory type loop