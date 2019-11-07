% fidCamera = fopen('/Users/tianqitang/CV_Research/images/meadow-3/dslr_calibration_undistorted/cameras.txt'); % TODO: Change directory
oneLine = fgetl(fidCamera);
counter = 1;
while(ischar(oneLine))
    if strcmp(oneLine(1), '#')
        oneLine = fgetl(fidCamera);
        continue;
    end
    space = find(oneLine == ' ');
    cameras.cameraId(counter) = str2num(oneLine(1:space(1)));
    cameras.camera{counter}.cameraModel = oneLine(space(1)+1:space(2));
    cameras.camera{counter}.imageWidth = oneLine(space(2)+1:space(3));
    cameras.camera{counter}.imageHeight = oneLine(space(3)+1:space(4));
    cameras.camera{counter}.fc(1) = str2num(oneLine(space(4)+1:space(5)));
    cameras.camera{counter}.fc(2) = str2num(oneLine(space(5)+1:space(6)));
    cameras.camera{counter}.cc(1) = str2num(oneLine(space(6)+1:space(7)));
    cameras.camera{counter}.cc(2) = str2num(oneLine(space(7)+1:end));
    counter = counter + 1;
    oneLine = fgetl(fidCamera);
end
fclose(fidCamera);

% fidImages = fopen('/Users/tianqitang/CV_Research/images/meadow-3/dslr_calibration_undistorted/images.txt'); % TODO: Change directory
oneLine = fgetl(fidImages);
counter = 1;
while(ischar(oneLine))
    if strcmp(oneLine(1), '#')
        oneLine = fgetl(fidImages);
        continue;
    end
    space = find(oneLine == ' ');
    image.id = str2num(oneLine(1:space(1)));
    image.quat(1) = str2num(oneLine(space(1)+1:space(2)));
    image.quat(2) = str2num(oneLine(space(2)+1:space(3)));
    image.quat(3) = str2num(oneLine(space(3)+1:space(4)));
    image.quat(4) = str2num(oneLine(space(4)+1:space(5)));
    image.T(1) = str2num(oneLine(space(5)+1:space(6)));
    image.T(2) = str2num(oneLine(space(6)+1:space(7)));
    image.T(3) = str2num(oneLine(space(7)+1:space(8)));
    image.camera = str2num(oneLine(space(8)+1:space(9)));
    image.path = oneLine(space(9)+1:end);
    images(counter) = image;
    oneLine = fgetl(fidImages);
    oneLine = fgetl(fidImages);
    counter = counter + 1;
end
fclose(fidImages);

% TODO: For each camera replace camera{1} as camera{id} image.camera will tell you which camera id was used to capture that image
cam1_num=find(cameras.cameraId==images(viewSelection(1)).camera);
cam2_num=find(cameras.cameraId==images(viewSelection(2)).camera);
K1 = [cameras.camera{cam1_num}.fc(1) 0 cameras.camera{cam1_num}.cc(1);...
     0 cameras.camera{cam1_num}.fc(2) cameras.camera{cam1_num}.cc(2);...
     0 0 1];
K2 = [cameras.camera{cam2_num}.fc(1) 0 cameras.camera{cam2_num}.cc(1);...
     0 cameras.camera{cam2_num}.fc(2) cameras.camera{cam2_num}.cc(2);...
     0 0 1];

% TODO: select different images by changing "viewSelection" you may need to manually select images that have enough overlap
R1 = quat2rotmSelf(images(viewSelection(1)).quat);
R2 = quat2rotmSelf(images(viewSelection(2)).quat);
T1 = transpose(images(viewSelection(1)).T);
T2 = transpose(images(viewSelection(2)).T);
C1 = -inv(R1) * T1;
C2 = -inv(R2) * T2;
R12 = R2 * R1';
T12 = R2 * (C1 - C2);

function R = quat2rotmSelf( q )
q2 = q';

s = q2(1,1);
x = q2(2,1);
y = q2(3,1);
z = q2(4,1);

% Explicitly define concatenation dimension for codegen
tempR = cat(1, 1 - 2*(y.^2 + z.^2),   2*(x.*y - s.*z),   2*(x.*z + s.*y),...
2*(x.*y + s.*z), 1 - 2*(x.^2 + z.^2),   2*(y.*z - s.*x),...
2*(x.*z - s.*y),   2*(y.*z + s.*x), 1 - 2*(x.^2 + y.^2) );

R = reshape(tempR, [3, 3, length(s)]);
R = permute(R, [2 1 3]);

end
