function videoMaker(frameMatrix,videoname,framerate)
    % VideoWriter オブジェクトを作成
    v = VideoWriter(strcat(videoname,'.mp4'), 'MPEG-4');
    % 時間区切りからフレームレートの計算と適用
    v.FrameRate = framerate;
    % 保存する動画の画質。数字の大きいほうが高画質.[0~100]
    v.Quality = 95;
    % ビデオの書き込みを開始
    open(v);
    sizedata = size(frameMatrix(1).cdata);
    for i = 1:numel(frameMatrix)
        immod = imresize(frameMatrix(i).cdata,sizedata(1:2));
        writeVideo(v,immod);
    end
    
    close(v);
end
