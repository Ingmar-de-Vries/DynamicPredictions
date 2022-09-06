function ActionPrediction_checkAtlases(cfg)

% check if correct atlases will be loaded for parcel definitions
addpath('\\XXX');
cfg.path = '\\XXX'; 

indir = fullfile(cfg.path,'data','MEG','PreprocessedSensor');
subfilz = dir(fullfile(indir,'preprocessedSensor*'));

for iSub = cfg.SubVec
    
    % source data
    indirSource = fullfile(cfg.path,'data','brainstorm_database','ActionPrediction','data',subfilz(iSub).name(end-7:end-4),[subfilz(iSub).name(end-7:end-4) '_allSensorData']);
    indirAtlas = fullfile(cfg.path,'data','brainstorm_database','ActionPrediction','anat',subfilz(iSub).name(end-7:end-4));
    
    % atlases
    sourcefile = dir(fullfile(indirSource,'*MEG_GRAD_MEG_MAG_KERNEL*'));% pick source files from Brainstorm database (source reconstruction was done in Brainstorm)
    fn2load = sprintf('%s%c%s',indirSource,filesep,sourcefile.name);
    kernel = load(fn2load);
    cortex = kernel.SurfaceFile(6:end-4);
    
    atlasfile = dir(fullfile(indirAtlas,'*pial_low*'));
    fn2load = sprintf('%s%c%s',indirAtlas,filesep,atlasfile(1).name);
    
    % check if the here loaded source file used the same atlas as the here loaded atlas file:
    if ~contains(fn2load,cortex)
        error(['cortical surface used for source reconstruction does not match the one in atlas for subject ' num2str(iSub)]);
    end
    
    atlas = load(fn2load);
    atlasIdx = logical(contains(extractfield(atlas.Atlas,'Name'),cfg.atlas) + contains(extractfield(atlas.Atlas,'Name'),'Schaefer_200'));
    atlas = atlas.Atlas(atlasIdx).Scouts;
    
    idx2remove = contains(extractfield(atlas,'Label'),'Background');
    atlas(idx2remove) = [];
    
    % fix subject 12 atlas such that it matches the one from all other subject
    % i.e., at some point during my analyses Brainstorm started using new Schaefer atlas with slightly different names and order for some of the parcels
    % I only received the individual anatomical scan for subject 12 after this had happened, so I had to change the parcels of subject 12 to the ones in all
    % other subjects, which is what is stored in the new2oldSchaeferXXX variable below
    if iSub == 12
        sub12 = extractfield(atlas,'Label')';
        allsubs = extractfield(referenceAtlas,'Label')';

        % first order according to hemisphere (left -- right) 
        lefties = false(size(sub12));
        for iroi = 1:length(atlas)
            if strcmp(sub12{iroi}(end),'L')
                lefties(iroi) = true;
            end
            allsubs{iroi}(1:14)=[];
        end
        atlas = [atlas(lefties) atlas(~lefties)];% and the actual atlas
        sub12 = [sub12(lefties) ; sub12(~lefties)];% temporary atlas labels
        
        % then manually change the mismatching names:
        new2oldSchaefer200 = [1:5 9 6:8 10:19 24:26 20:23 27:42 44:48 43 49:53 58 59 54:57 60 61 64 62 63 65:67 69 70 68 71:80 84 85 81 86 82 ...
            83 87:90 94 91:93 95:97 100 98 99 101:110 112:114 111 115:120 125 121:124 126:132 134:138 133 139:143 148:151 144:147 152 157 158 ...
            153 155 156 154 159:161 165 164 166 162 163 167:179 181:184 180 185:190 194 191:193 195 196 200 197:199];
        save('\\XXX\new2oldSchaefer200','new2oldSchaefer200');
        
        new2oldSchaefer400 = [1:30 40:46 31:39 47:100 108:112 101:107 113:115 118:121 116 117 122:127 131:135 128:130 136:155 166 167 159:163 ...
            156:158 164 165 168:170 171:177 187 178 179 180:184 188 185 186 189 196 197 198 190:193 199 200 194 195 201:235 245:249 236:244 ...
            250:266 269:280 267 268 281:292 299:304 293:298 312:315 309:311 305 316 317 306 318 319 307 308 320:324 330 327:329 331 332 325 ...
            326 333:353 365:367 356:358 354 359:361 355 362:364 368:378 388 379 380:385 389 386 387 390 397 398 391:394 399 400 395 396];
        save('\\XXX\new2oldSchaefer400','new2oldSchaefer400');
    end
    
    % set atlas of subject 1 as reference atlas and compare all others
    % against it to make sure we're using the same atlas for all subjects
    if iSub == 1
        referenceAtlas = atlas;
    else
        if ~all(strcmp(extractfield(referenceAtlas,'Label'),extractfield(atlas,'Label')))
            error(['The atlas of subject ' num2str(iSub) ' does not match the reference atlas'])
        end
    end
    
end


end