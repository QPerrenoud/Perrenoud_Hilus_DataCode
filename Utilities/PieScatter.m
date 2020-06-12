function PieScatter(db1X, db1Y, db2P, db2Color, db1XLim, db1YLim)
narginchk(3, 6)
% 
% % Creates bogus input arguments for testing
% db1X = rand(10, 1);
% db1Y = rand(10, 1);
% db2P = round(rand(10, 3));
% nargin = 3;
% figure('Position', [100 100 1600 900])

% Checks that the compulsory input arguments are properly formatted
if ~isvector(db1X) || ~isvector(db1Y); error('db1X and db1Y must be vectors'); end
if isrow(db1X), db1X = db1X'; end
if isrow(db1Y), db1Y = db1Y'; end
if length(db1X) ~= length(db1Y); error('db1X and db1Y must have equal length'); end
if size(db2P, 1) ~= length(db1Y); error('dbXP must have as many rows as element in db1X and db1Y'); end

% Initializes optional input arguments if they don't exist
if nargin < 4; db2Color = []; end
if nargin < 5; db1XLim = []; end
if nargin < 6; db1YLim = []; end
    
% Check that optional input arguments are properly formatted.
if any(size(db2Color) ~= [size(db2P, 2), 3])
    db2Jet      = jet;
    db2Color    = db2Jet(round(linspace(1, 64, size(db2P, 2))), :); 
end
if any(size(db1XLim) ~= [1 2])
    db1XLim = [min(db1X) - .1 * range(db1X) max(db1X) + .1 * range(db1X)];
end
if any(size(db1YLim) ~= [1 2])
    db1YLim = [min(db1Y) - .1 * range(db1Y) max(db1Y) + .1 * range(db1Y)];
end

% Calculate the size of the markers in x an y
    % gets the aspect ratio and the long axis
db1FPos = get(gcf, 'Position');
db1APos = get(gca, 'Position');
db1AR = (db1APos(3) * db1FPos(3)) ./ (db1APos(4) * db1FPos(4));
    % calculates the size in x and y
dbSz_X = range(db1XLim) * 0.01;
dbSz_Y = range(db1YLim) * 0.01 * db1AR;

% Computes the marker coordinates relative to the point of interest
db1XM = cosd(0:360) * dbSz_X; 
db1YM = sind(0:360) * dbSz_Y;

%Calculates useful dimensions
inNPt   = length(db1X);
inNPie  = size(db2P, 2);

% Calculates the indices of pies out of proportion for each point
db2PieBnd = [sum(db2P, 2) ~= 0 round((cumsum(db2P, 2)./sum(db2P, 2)) * 360) + 1];

% Loops through points
hold on
for iPt = 1:inNPt
    if db2PieBnd(iPt, 1) == 0
        fill(db1XM + db1X(iPt), db1YM + db1Y(iPt), [0 0 0], 'LineStyle', 'none')
    else
        for iCnd = 1:inNPie
            in1Idx = db2PieBnd(iPt, iCnd):db2PieBnd(iPt, iCnd + 1);
            db1XM_Pie = [0 db1XM(in1Idx)] + db1X(iPt);
            db1YM_Pie = [0 db1YM(in1Idx)] + db1Y(iPt);
            fill(db1XM_Pie, db1YM_Pie, db2Color(iCnd, :), 'LineStyle', 'none')
        end
    end
end
xlim(db1XLim), ylim(db1YLim);