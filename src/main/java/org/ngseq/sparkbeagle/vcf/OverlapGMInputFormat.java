package org.ngseq.sparkbeagle.vcf;

/**
 * Created by maarala1 on 11/7/18.
 */

import htsjdk.tribble.readers.TabixReader;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.*;
import org.apache.hadoop.mapreduce.lib.input.FileSplit;
import org.apache.hadoop.mapreduce.lib.input.LineRecordReader;
import org.apache.hadoop.mapreduce.lib.input.TextInputFormat;
import org.apache.hadoop.util.LineReader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class OverlapGMInputFormat extends TextInputFormat {
    public static final String LINES_PER_MAP = "mapreduce.input.lineinputformat.linespermap";

    private static Path filePath;
    private static TabixReader tabixReader;
    private static List<String> intervals = new ArrayList<>();

    public static void setTabixReader(TabixReader tabixReader) {
        OverlapGMInputFormat.tabixReader = tabixReader;
    }

    public static void setIntervals(List<String> intervals) {
        OverlapGMInputFormat.intervals = intervals;
    }

    public static void setFilePath(Path path) {
        OverlapGMInputFormat.filePath = path;
    }

    public OverlapGMInputFormat() {
    }

    public RecordReader<LongWritable, Text> createRecordReader(InputSplit genericSplit, TaskAttemptContext context) {
        context.setStatus(genericSplit.toString());
        return new LineRecordReader();
    }

    @Override
    public List<InputSplit> getSplits(JobContext job) throws IOException {
        List<InputSplit> splits = new ArrayList();
        int numLinesPerSplit = getNumLinesPerSplit(job);
        int overlap = getOverlapPerSplit(job);
        int offset = getOffset(job);
        Iterator it = this.listStatus(job).iterator();

        //TODO:get positions from Genmap file




        while(it.hasNext()) {
            FileStatus status = (FileStatus) it.next();
            splits.addAll(getSplitsForFile(status, job.getConfiguration(), OverlapGMInputFormat.intervals, overlap, offset));
        }

        return splits;
    }

    public static List<FileSplit> getSplitsForFile(FileStatus status, Configuration conf, List<String> intervals, int overlap, int offset) throws IOException {
        List<FileSplit> splits = new ArrayList();
        Path fileName = status.getPath();
        long len = status.getLen();
        if(status.isDirectory()) {
            throw new IOException("Not a file: " + fileName);
        } else {
            FileSystem fs = fileName.getFileSystem(conf);
            LineReader lr = null;

            try {
                //TODO: uncompress here if needed
                FSDataInputStream in = fs.open(fileName);
                lr = new LineReader(in, conf);
                Text line = new Text();
                int numLines = 0;
                long begin = 0L;
                long length = 0L;
                long lensofar = 0L;
                long olLen = 0L;

                int num;
                boolean isfirst = true;
                int splitcnt = 0;
                while((num = lr.readLine(line)) > 0) {
                    ++numLines;

                    lensofar += (long)num;

                    //TODO:determine overlap length some other way
                    /*int rest = numLinesPerSplit-numLines;
                    if(rest < overlap) {
                        olLen += (long)num;
                    }*/
                    length += (long)num;

                    if(line.toString().split("\t")[1].equals(intervals.get(splitcnt))) {

                        System.out.println("length:"+length);
                        System.out.println("lensofar:"+lensofar);
                        System.out.println("overlap:"+olLen);
                        System.out.println("begin:"+begin);
                        System.out.println("numlines:"+numLines);
                        System.out.println("filename:"+fileName.toString());

                        if(!isfirst)
                            splits.add(createFileSplit(fileName, begin, lensofar-begin));
                        else
                            splits.add(createFileSplit(fileName, 0L, lensofar));

                        begin = (lensofar-olLen);

                        olLen = 0L;
                        length = 0L;
                        numLines = 0;
                        isfirst=false;
                        splitcnt++;
                    }

                    /*if(numLines == numLinesPerSplit+overlap ) {

                    }*/

                }

                if(numLines != 0) {

                    System.out.println("_file len:"+len);
                    System.out.println("_begin+length:"+(begin+length));
                    System.out.println("_begin:"+begin);
                    System.out.println("_lensofar:"+lensofar);
                    System.out.println("_length:"+length);
                    System.out.println("_numlines_:"+numLines);
                    System.out.println("_remaining len:"+(len-begin));
                    System.out.println("filename:"+fileName.toString());

                    splits.add(createFileSplit(fileName, begin, lensofar-begin));

                }
            } finally {
                if(lr != null) {
                    lr.close();
                }

            }

            return splits;
        }
    }

    protected static FileSplit createFileSplit(Path fileName, long begin, long length) {
        return begin == 0L?new FileSplit(fileName, begin, length - 1L, new String[0]):new FileSplit(fileName, begin - 1L, length, new String[0]);
    }

    public static void setNumLinesPerSplit(Job job, int numLines) {
        job.getConfiguration().setInt("mapreduce.input.lineinputformat.linespermap", numLines);
    }
    public static void setOverlapPerSplit(Job job, int overlap) {
        job.getConfiguration().setInt("mapreduce.input.lineinputformat.overlaplines", overlap);
    }
    public static void setOffset(Job job, int offset) {
        job.getConfiguration().setInt("mapreduce.input.lineinputformat.offset", offset);
    }
    public static int getOverlapPerSplit(JobContext job) {
        return job.getConfiguration().getInt("mapreduce.input.lineinputformat.overlaplines", 1);
    }

    public static int getOffset(JobContext job) {
        return job.getConfiguration().getInt("mapreduce.input.lineinputformat.offset", 0);
    }

    public static int getNumLinesPerSplit(JobContext job) {
        return job.getConfiguration().getInt("mapreduce.input.lineinputformat.linespermap", 1);
    }
}

