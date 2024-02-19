for cell_num = 3
    i = 0;
    outputPath = 'Output/';
    out = []; params = [];
for j = 1:10
    Dataset = import_parquet(sprintf('Cell%d_RPT%d.parquet', cell_num, i));
    cycle= 0;
    reset_path;
    [out, params] = main_one_multi(Dataset, out, params, cycle, outputPath, j);
end
end