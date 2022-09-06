function DynamicPredictions_defineSourceROIs()

%%  define ROIs as parcels or combinations of parcels from the HCP atlas

% visual areas
ROIdefinition.parcels{1} = {'V1'};
ROIdefinition.names{1} = {'V1'}; 

ROIdefinition.parcels{2} = {'V2'};
ROIdefinition.names{2} = {'V2'};

ROIdefinition.parcels{3} = {'V3','V4'};
ROIdefinition.names{3} = {'V3V4'};% combined to keep #vertices similar across ROIs and because MEG can't clearly distinguish between V3 and V4 anyway 

% action observation network (AON)
ROIdefinition.parcels{4} = {'V4t', 'FST', 'MT', 'MST', 'LO1', 'LO2', 'LO3', 'PH', 'PHT', 'TPOJ2', 'TPOJ3'};
ROIdefinition.names{4} = {'LOTC'};%posterior Lateral Occipital Temporal Cortex

ROIdefinition.parcels{5} = {'PF', 'PFt', 'AIP', 'IP2'};
ROIdefinition.names{5} = {'aIPL'};% anterior Inferior Parietal Lobule

ROIdefinition.parcels{6} = {'IFJa', 'IFJp', '6r', '6v', 'PEF', 'IFSp', '44', '45'};
ROIdefinition.names{6} = {'VentPremot'};%

save('\\XXX\ROIdefinitions','ROIdefinition');
