
%%
for phase = 0 : 4
    SLM = patterntabArb([2, -11], 4.064, 0.334, phase, 0, 5, 1920, 1152);
    SLM = SLM + 21;

    imwrite(SLM, strcat('image0', num2str(phase), '.bmp'),'BMP');
end

%%
for phase = 5 : 9
    SLM = patterntabArb([14, -5], 4.075, 0.338, phase, 0, 5, 1920, 1152);
    SLM = SLM + 50;

    imwrite(SLM, strcat('image0', num2str(phase), '.bmp'),'BMP');
end

%%
for phase = 10 : 14
    SLM = patterntabArb([-13, -11], 4.067, 0.336, phase, 0, 5, 1920, 1152);
    SLM = SLM + 44;

    imwrite(SLM, strcat('image', num2str(phase), '.bmp'),'BMP');
end


%%
function SLM = patterntabArb(vecA, period, onfrac, phaseInd, phaseOffset, nphases, sizex, sizey)

%{
1. Decide on a vector parallel to the pattern stripes (vecA) and the period
2. Obtain the vector perpendicular to vecA that has the length of the desired period (vecB)
3. For any pixel (x, y), calculate the projection (or inner product) of vector {x, y} onto vecB, take its modulo (remainder after dividing by the period) value because of the periodicity, and judge if it’s greater or less than the duty cycle (onfrac) and assign 0 and 1 accordingly.
4. Repeat the above for the other phase steps simply by moving the origin of vecB along vecB by the length of the phase step.
%}

veckA = [1, -1] .* flip(vecA);                          % veckA = [-15, -4]
vecB = veckA ./ sqrt(vecA(1)^2 + vecA(2)^2) * period;   % vecB = [-4.8312, -1.2883]
area = vecB(1) * veckA(1) + vecB(2) * veckA(2);         % area = sqrt(15 ^ 2 + 4 ^ 2) * 5 = 77.6209
onpix = area * onfrac;                                  % onpix = 23.2863
phaseStep = vecB / nphases;                             % phaseStep = [-0.9662, -0.2577]

SLM = zeros([sizey, sizex], 'uint8');

for i = 0 : (sizex - 1)
    for j = 0 : (sizey - 1)
        pixel_ind = [i, j] - (phaseStep * phaseInd + phaseOffset /(2 * pi) * vecB); 
        pixel_vecB = pixel_ind(1) * veckA(1) + pixel_ind(2) * veckA(2);
        
        %{
        Suppose phaseInd = 0, phaseOffset = 0:
        (i, j) = (0, 0), pixel_ind = [0, 0], pixel_vecB = 0, mod(0, 77.6209) = 0 < 23.2863, SLM(1, 1) = 128
        (i, j) = (3, 3), pixel_ind = [3, 3], pixel_vecB = -57, mod(-57, 77.6209) = 20.6209 < 23.2863, SLM(4, 4) = 128
        (i, j) = (5, 5), pixel_ind = [5, 5], pixel_vecB = -95, mod(-95, 77.6209) = 60.2418 >= 23.2863, SLM(4, 4) = 0
        
        Suppose phaseInd = 1, phaseOffset = 0:
        (i, j) = (0, 0), pixel_ind = [0.9662, 0.2577], pixel_vecB = -15.5242, mod(-15.5242, 77.6209) = 62.0967 >= 23.2863, SLM(1, 1) = 0
        (i, j) = (3, 3), pixel_ind = [3.9662, 3.2577], pixel_vecB = -72.5242, mod(-72.5242, 77.6209) = 5.0967 < 23.2863, SLM(4, 4) = 128
        (i, j) = (5, 5), pixel_ind = [5.9662, 5.9662], pixel_vecB = -110.5242, mod(-110.5242, 77.6209) =  44.7176 >= 23.2863, SLM(4, 4) = 0
        %}
        
        if mod(pixel_vecB, area) >= onpix
            SLM(j + 1, i + 1) = 0;
        else
            SLM(j + 1, i + 1) = 128;
        end
    end
end

end