function genomeObject = FastaToGenomeObject(regionFasta,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'maxProbesPerRegion','integer',inf};
defaults(end+1,:) = {'probeLength','integer',40};
parameters = ParseVariableArguments(varargin, defaults, mfilename);

regionFastaData = fastaread(regionFasta);

% hack the Transcriptome object that requires gene='' in the name
fastaHeader = {regionFastaData.Header}';
nEntries = length(fastaHeader);
if ~contains(fastaHeader{1},'gene=')
    for i=1:nEntries
        currName = regionFastaData(i).Header;
        regionFastaData(i).Header = ['',' gene=',currName];
        maxLength = parameters.probeLength*parameters.maxProbesPerRegion*2;
        if length(regionFastaData(i).Sequence) > maxLength
           regionFastaData(i).Sequence = regionFastaData(i).Sequence(1:maxLength); 
        end
    end
end

genomeObject = Transcriptome(regionFastaData,'verbose',false);