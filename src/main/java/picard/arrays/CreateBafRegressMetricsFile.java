package picard.arrays;

/*
 * The MIT License
 *
 * Copyright (c) 2020 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A simple program to create a standard picard metrics file
 * from the output of bafRegress
 */
@CommandLineProgramProperties(
        summary = CreateBafRegressMetricsFile.USAGE_DETAILS,
        oneLineSummary = "Program to generate a picard metrics file from the output of the bafRegress tool.",
        programGroup = picard.cmdline.programgroups.GenotypingArraysProgramGroup.class
)

@DocumentedFeature
public class CreateBafRegressMetricsFile extends CommandLineProgram {
    static final String USAGE_DETAILS =
            "CreateBafRegressMetricsFile takes an output file as generated by the bafRegress tool and creates a picard metrics file. " +
                    "BAFRegress <a href='https://genome.sph.umich.edu/wiki/BAFRegress'>bafRegress</a> " +
                    "is a software that detects and estimates sample contamination using B allele frequency data from Illumina genotyping arrays using a regression model." +
                    "<h4>Usage example:</h4>" +
                    "<pre>" +
                    "java -jar picard.jar CreateBafRegressMetricsFile \\<br />" +
                    "      INPUT=bafRegress.output.txt \\<br />" +
                    "      OUTPUT=outputBaseName" +
                    "</pre>";

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The output of bafRegress (typically captured stdout).")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Basename for the metrics file that will be written." +
            " Resulting file will be <OUTPUT>." + FILE_EXTENSION)
    public File OUTPUT;

    public static final String FILE_EXTENSION = "bafregress_metrics";

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        final File metricsFile = new File(OUTPUT + "." + FILE_EXTENSION);
        IOUtil.assertFileIsWritable(metricsFile);

        final MetricsFile<BafRegressMetrics, ?> bafRegressMetricsMetricsFile = getMetricsFile();

        try (BufferedReader br = new BufferedReader(new FileReader(INPUT))) {
            String line;
            line = br.readLine();
            if (!line.equals("sample\testimate\tstderr\ttval\tpval\tcallrate\tNhom")) {
                throw new PicardException("Unrecognized header line: '" + line + "' in " + INPUT.getAbsolutePath());
            }
            while ((line = br.readLine()) != null) {
                String[] entries = line.split("\\s+");
                if (entries.length != 7) {
                    throw new IOException("Invalid number of entries (" + entries.length + ") in line: " + line);
                }
                // Load up and store the metrics
                final BafRegressMetrics metrics = new BafRegressMetrics();
                metrics.SAMPLE = entries[0];
                metrics.ESTIMATE = Double.parseDouble(entries[1]);
                metrics.STDERR = Double.parseDouble(entries[2]);
                metrics.TVAL = Double.parseDouble(entries[3]);
                metrics.PVAL = Double.parseDouble(entries[4]);
                metrics.LOG10_PVAL = Math.log10(metrics.PVAL);
                metrics.CALL_RATE = Double.parseDouble(entries[5]);
                metrics.NHOM = Integer.parseInt(entries[6]);

                bafRegressMetricsMetricsFile.addMetric(metrics);
            }
            bafRegressMetricsMetricsFile.write(metricsFile);
        } catch (IOException e) {
            throw new PicardException("Error parsing bafRegress Output", e);
        }
        return 0;
    }
}

