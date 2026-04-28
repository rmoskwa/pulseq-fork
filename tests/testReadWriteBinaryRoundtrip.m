function tests = testReadWriteBinaryRoundtrip
    tests = functiontests(localfunctions);
end

function test_seq_text_binary_text_roundtrip(testCase)
    currDir = fileparts(mfilename('fullpath'));

    try
        mr.opts();
    catch
        addpath(fullfile(currDir,'../matlab'));
    end

    seqDir = getenv('PULSEQ_SEQ_DIR');
    if isempty(seqDir)
        %seqDir = fullfile(currDir, '../matlab');
        seqDir = currDir;
    end

    listing = dir(fullfile(seqDir, '*.seq'));
    testCase.assertNotEmpty(listing, sprintf('No .seq files found in directory: %s', seqDir));

    tmpRoot = tempname;
    mkdir(tmpRoot);
    cleanupObj = onCleanup(@() localCleanup(tmpRoot)); %#ok<NASGU>

    for k = 1:numel(listing)
        srcPath = fullfile(listing(k).folder, listing(k).name);
        fprintf('Roundtrip check: %s\n', listing(k).name);
        [~, baseName, ~] = fileparts(srcPath);
        binPath = fullfile(tmpRoot, [baseName '.bin']);
        canonicalSeqPath = fullfile(tmpRoot, [baseName '_canonical.seq']);
        outSeqPath = fullfile(tmpRoot, [baseName '_roundtrip.seq']);

        seq = mr.Sequence();
        seq.read(srcPath);
        seq.write(canonicalSeqPath);
        seq.writeBinary(binPath);

        seqRound = mr.Sequence();
        seqRound.readBinary(binPath);
        seqRound.write(outSeqPath);

        srcTxt = readAndNormalizeSeq(srcPath);
        canonicalTxt = readAndNormalizeSeq(canonicalSeqPath);
        outTxt = readAndNormalizeSeq(outSeqPath);
        if ~strcmp(outTxt, srcTxt)
            % Some reference files are not in canonical formatting and may
            % change under read()->write() even without binary conversion.
            % In that case, compare against canonical text generated from
            % the same loaded sequence object.
            if strcmp(outTxt, canonicalTxt)
                continue;
            end
            roundtrip_mismatch_source=fullfile(tmpRoot,'roundtrip_mismatch_source.seq');
            roundtrip_mismatch_canonical=fullfile(tmpRoot,'roundtrip_mismatch_canonical.seq');
            roundtrip_mismatch_output=fullfile(tmpRoot,'roundtrip_mismatch_output.seq');
            writeTextFile(roundtrip_mismatch_source, srcTxt);
            writeTextFile(roundtrip_mismatch_canonical, canonicalTxt);
            writeTextFile(roundtrip_mismatch_output, outTxt);
            diffPos = firstDiffPosition(srcTxt, outTxt);
            testCase.verifyFail(sprintf(['Roundtrip mismatch for file: %s. ' ...
                'First differing position: %d. Debug copies written to: ' ...
                roundtrip_mismatch_source ', ' roundtrip_mismatch_canonical ...
                ' and ' roundtrip_mismatch_output], ...
                listing(k).name, diffPos));
        end
    end
end

function txt = readAndNormalizeSeq(path)
    txt = fileread(path);
    txt = strrep(txt, sprintf('\r\n'), sprintf('\n'));
    txt = strrep(txt, sprintf('\r'), sprintf('\n'));
    txt = regexprep(txt, '\n\[SIGNATURE\][\s\S]*$', '');
    txt = regexprep(txt, '[ \t]+\n', sprintf('\n'));
    txt = regexprep(txt, '\n+$', sprintf('\n'));
end

function localCleanup(path)
    if exist(path, 'dir')
        try
            rmdir(path, 's');
        catch
        end
    end
end

function writeTextFile(path, txt)
    fid = fopen(path, 'w');
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fwrite(fid, txt);
end

function pos = firstDiffPosition(a, b)
    n = min(length(a), length(b));
    idx = find(a(1:n) ~= b(1:n), 1, 'first');
    if isempty(idx)
        pos = n + 1;
    else
        pos = idx;
    end
end
