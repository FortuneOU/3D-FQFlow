%% =======================
%  参数设置
% ========================
resultFolder = './Data/';
matFiles = dir(fullfile(resultFolder,'CmodeImages*.mat')); % 根据需要调整
outputVideo = fullfile(resultFolder,'DopplerVideo.mp4');



% 图像坐标
x_image = squeeze(xi(1,:,1)) * 100;  % cm
y_image = squeeze(zi(1,1,:)) * 100;  % cm

% =======================
%  创建视频对象
% =======================
v = VideoWriter('DopplerVideo.mp4','MPEG-4');
v.FrameRate = 10;
open(v);

% =======================
%  循环读取每个 mat 中的 Vd，并保存为视频
% =======================
figure;
for f = 1:length(matFiles)

    % 加载当前 mat 文件
    data = load(fullfile(matFiles(f).folder, matFiles(f).name));
    Vd = data.Vd;    % 维度应为 [x, z, frame]

    numFrames = size(Vd,3);

    for k = 1:numFrames
        imagesc(x_image, y_image, Vd(:,:,k));
        colormap(dopplermap);
        colorbar;
        caxis([-0.6 0.6] * max(abs(Vd(:))));
        axis equal tight;
        title(sprintf('Frame %d of file %d',k,f));

        frame = getframe(gcf);
        writeVideo(v, frame);
    end
end

close(v);
disp('视频已保存完成');


figure;imagesc( Vd(:,:,40))
