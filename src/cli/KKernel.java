package cli;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.converters.FileConverter;

import benchmark.dataset.AttrDataset;
import benchmark.dataset.SDataset;
import cli.KCommon.RangeConverter;
import cli.KCommon.RangeSplitter;
import comparison.kernel.ExplicitMappingKernel;
import comparison.kernel.Kernel;
import comparison.kernel.basic.DiracKernel;
import comparison.kernel.graph.EdgeKernel;
import comparison.kernel.graph.FixedLengthRandomWalkKernel;
import comparison.kernel.graph.NodeCentricWalkKernel;
import comparison.kernel.graph.GraphHopper;
import comparison.kernel.graph.GraphletKernel;
import comparison.kernel.graph.MaxFixedLengthRandomWalkKernel;
import comparison.kernel.graph.VertexKernel;
import comparison.kernel.graph.WeisfeilerLehmanSubtreeKernel;
import graph.LGraph;
import util.kernel.GramUtil;


public class KKernel {
	
	static class CommandMain extends KCommon.CommandMain {

		@Parameter(names = { "-D", "--datadir" }, description = "Directory containing the data files", converter = FileConverter.class)
		File dataDir = new File("data");
		
		@Parameter(names = { "-G", "--gramdir" }, description = "Output directory where the gram files are stored", converter = FileConverter.class)
		File gramDir = new File("gram");

		@Parameter(names = { "-e", "--explicit" }, description = "Force computation by explicit feature maps")
		boolean explicit = false;

		@Parameter(names = { "-i", "--implicit" }, description = "Force computation by implicit feature maps")
		boolean implicit = false;
		
		@Parameter(names = { "-l", "--log" }, description = "File for saving running times", converter = FileConverter.class)
		File logFile = new File("log_runtime.txt");
		
		@Parameter(names = {"-d", "--datasets"}, description = "List of data sets")
		List<String> datasets;
		
		@Parameter(names = { "-a", "--all" }, description = "Compute kernel for all data sets in the data directory")
		private boolean all;

	}
	
	public static abstract class KernelConfig {
		
		abstract SDataset preprocessDataset(AttrDataset ds);
		
		abstract ArrayList<Kernel<LGraph<String, String>>> getKernels();
	}
	
	public static abstract class SimpleGraphKernelConfig extends KernelConfig {
		
		SDataset preprocessDataset(AttrDataset ds) {
			return ds.getSDataset();
		}
	}
	
	public static abstract class CommandWL extends SimpleGraphKernelConfig {

		@Parameter(names = { "-h", "--height" }, 
				description = "height, i.e., the number of refinement steps", 
				variableArity = true,
				splitter = RangeSplitter.class,
				converter = RangeConverter.class)
		List<Integer> iterations = Arrays.asList(0,1,2,3,4,5);

		ArrayList<Kernel<LGraph<String, String>>> getKernels() {
			ArrayList<Kernel<LGraph<String, String>>> kernels = new ArrayList<>();
			for (Number h : iterations) {
				kernels.add(getKernel(h.intValue()));
			}
			return kernels;
		}

		abstract Kernel<LGraph<String, String>> getKernel(int height);
	}
	
	@Parameters(commandDescription = "Compute the Weisfeiler-Lehman subtree kernel.")
	public static class CommandWLS extends CommandWL {
		Kernel<LGraph<String, String>> getKernel(int height) {
			return new WeisfeilerLehmanSubtreeKernel<>(height);
		}
	}
	
	@Parameters(commandDescription = "Compute the GraphHopper kernel, where vertex labels are compared by the "
			+ "Dirac kernel.")
	public static class CommandGH extends SimpleGraphKernelConfig {
		@Override
		ArrayList<Kernel<LGraph<String, String>>> getKernels() {
			ArrayList<Kernel<LGraph<String, String>>> kernels = new ArrayList<>();
			DiracKernel dk = new DiracKernel();
			kernels.add(new GraphHopper<String,String>(dk, false));
			return kernels;
		}
	}

	public abstract static class CommandFLRWAbstact extends SimpleGraphKernelConfig {
		
		@Parameter(names = { "-l", "--length" }, 
				description = "walk length", 
				variableArity = true,
				splitter = RangeSplitter.class,
				converter = RangeConverter.class)
		List<Integer> lengths = Arrays.asList(0, 1,2,3,4,5);
		
		@Override
		ArrayList<Kernel<LGraph<String, String>>> getKernels() {
			ArrayList<Kernel<LGraph<String, String>>> kernels = new ArrayList<>();
			for (int length: lengths) {
				kernels.add(getKernel(length));
			}
			return kernels;
		}
		
		abstract Kernel<LGraph<String, String>> getKernel(int length);
	}
	
	@Parameters(commandDescription = "Compute the fixed length random walk kernel, where "
			+ "vertex and edge labels are compared by the Dirac kernel.")
	public static class CommandFLRW extends CommandFLRWAbstact {
		
		Kernel<LGraph<String, String>> getKernel(int length) {
			DiracKernel dk = new DiracKernel();
			return new FixedLengthRandomWalkKernel<String, String>(length, dk, dk);
		}		
	}
	
	@Parameters(commandDescription = "Compute the sum of random walk kernels up to a fixed length, "
			+ "where vertex and edge labels are compared by the Dirac kernel.")
	public static class CommandMFLRW extends CommandFLRWAbstact {
		
		Kernel<LGraph<String, String>> getKernel(int length) {
			DiracKernel dk = new DiracKernel();
			return new MaxFixedLengthRandomWalkKernel<String, String>(length, dk, dk);
		}		
	}
	
	@Parameters(commandDescription = "Compute the node-centric walk kernel, "
			+ "where vertex and edge labels are compared by the Dirac kernel.")
	public static class CommandNCW extends CommandFLRWAbstact {
		boolean wl = false;
		@Parameter(names = { "-a", "--alpha" }, 
				description = "strictness of neighborhood comparison; parameter of Gaussian RBF", 
				variableArity = true,
				splitter = RangeSplitter.class,
				converter = RangeConverter.class)
		List<Double> alphas = Arrays.asList(0.01, 0.1, 1.0, 10.0);
		@Parameter(names = { "-b", "--beta" }, 
				description = "weight of graph-level walk counts", 
				variableArity = true,
				splitter = RangeSplitter.class,
				converter = RangeConverter.class)
		List<Double> betas = Arrays.asList(0.0, 0.5, 1.0);

		
		ArrayList<Kernel<LGraph<String, String>>> getKernels() {
			DiracKernel dk = new DiracKernel();
			ArrayList<Kernel<LGraph<String, String>>> kernels = new ArrayList<>();
			for (int length: lengths) {
				for (double alpha : alphas) {
					for (double beta : betas) {
						kernels.add(new NodeCentricWalkKernel<String, String>(length, alpha, beta, dk, dk, wl));
					}
				}
			}
			return kernels;
		}

		@Override
		Kernel<LGraph<String, String>> getKernel(int length) {	return null; } // not uses, since two parameters must be considered		
	}
	
	@Parameters(commandDescription = "Compute the node-centric walk kernel with WL expressivity,"
			+ "where vertex and edge labels are compared by the Dirac kernel.")
	public static class CommandNCWWL extends CommandNCW {
		CommandNCWWL() { this.wl = true; }
	}
	
	@Parameters(commandDescription = "Compute the graphlet kernel taking connected induced "
			+ "subgraphs on three vertices and discrete vertex and edge labels into account.")
	public static class CommandGL3 extends SimpleGraphKernelConfig {

		@Override
		ArrayList<Kernel<LGraph<String, String>>> getKernels() {
			ArrayList<Kernel<LGraph<String, String>>> kernels = new ArrayList<>();
			kernels.add(new GraphletKernel<>());
			return kernels;
		}
	}

	@Parameters(commandDescription = "Compute the vertex label kernel.")
	public static class CommandVL extends SimpleGraphKernelConfig {

		@Override
		ArrayList<Kernel<LGraph<String, String>>> getKernels() {
			ArrayList<Kernel<LGraph<String, String>>> kernels = new ArrayList<>();
			kernels.add(new VertexKernel<>(new DiracKernel()));
			return kernels;
		}
	}
	
	@Parameters(commandDescription = "Compute the edge label kernel, where an edge label is "
			+ "considered a triple of the labels of the edge itself and its endpoints."
			+ "Each edge generates two triples with interchanged endpoints.")
	public static class CommandEL extends SimpleGraphKernelConfig {

		@Override
		ArrayList<Kernel<LGraph<String, String>>> getKernels() {
			ArrayList<Kernel<LGraph<String, String>>> kernels = new ArrayList<>();
			DiracKernel dk = new DiracKernel();
			kernels.add(new EdgeKernel<>(dk, dk));
			return kernels;
		}
	}
	
	static CommandMain cm = new CommandMain();
	static CommandWLS wls = new CommandWLS();
	static CommandGH gh = new CommandGH();
	static CommandFLRW flrw = new CommandFLRW();
	static CommandMFLRW mflrw = new CommandMFLRW();
	static CommandNCW ncw = new CommandNCW();
	static CommandNCWWL ncwwl = new CommandNCWWL();
	static CommandGL3 gl3 = new CommandGL3();
	static CommandVL vl = new CommandVL();
	static CommandEL el = new CommandEL();

	@SuppressWarnings({ "unchecked", "rawtypes" })
	public static void main(String[] argsString) throws IOException, InterruptedException {
		
		JCommander jc = new JCommander(cm);
		jc.addCommand("wls", wls);
		jc.addCommand("gh", gh);
		jc.addCommand("flrw", flrw);
		jc.addCommand("mflrw", mflrw);
		jc.addCommand("ncw", ncw);
		jc.addCommand("ncwwl", ncwwl);		
		jc.addCommand("gl3", gl3);
		jc.addCommand("vl", vl);
		jc.addCommand("el", el);
		
		jc.getMainParameter();
		jc.setProgramName("kkernel");

		jc.parse(argsString);

		if (cm.help || jc.getParsedCommand() == null) {
			jc.usage();
			System.exit(0);
		}
		
	    if (cm.explicit && cm.implicit) {
	        throw new ParameterException("Choose either -i or -e.");
	    }
	    
	    if (!cm.all && cm.datasets == null) {
	        throw new ParameterException("No datasets specified.");
	    }
		
		KernelConfig kc = null;

		switch (jc.getParsedCommand()) {
			case "wls"   : kc = wls;    break;
			case "gh"    : kc = gh;     break;
			case "flrw"  : kc = flrw;   break;
			case "mflrw" : kc = mflrw;  break;
			case "ncw"   : kc = ncw;    break;
			case "ncwwl" : kc = ncwwl;  break;
			case "gl3"   : kc = gl3;    break;
			case "vl"    : kc = vl;     break;
			case "el"    : kc = el;     break;
                }

		if (!cm.gramDir.exists()) cm.gramDir.mkdir();
		
		if (cm.all) {
			cm.datasets = KCommon.getDatasets(cm.dataDir);
		}
		
		for (String dName : cm.datasets) {
			AttrDataset ds = KCommon.load(dName, cm.dataDir);
			SDataset ds2 = kc.preprocessDataset(ds);
			for (Kernel<LGraph<String, String>> k : kc.getKernels()) {

				System.out.println("Kernel:   "+k.getID());
				System.out.println("Dataset:  "+ds2.getID()+"  converted from: "+ds.getID());
				double[][] gram;
				boolean explicit = false;
				long startTime = System.nanoTime();
				if (cm.explicit) {
					gram = ((ExplicitMappingKernel)k).computeExplicit(ds2);
					explicit = true;
				} else if (cm.implicit || !(k instanceof ExplicitMappingKernel)) {
					gram = k.computeAll(ds2);
					explicit = false;
				} else {
					ExplicitMappingKernel emk = (ExplicitMappingKernel)k;
					try {
						gram = emk.computeExplicit(ds2);
						explicit = true;
					} catch (IllegalStateException e) {
						System.out.println("Non-explicit computation due to kernel choice!");
						gram = k.computeAll(ds2);
						explicit = false;
					}
				}
				long runtime = System.nanoTime() - startTime;
				
				// write running time
				FileWriter fw = new FileWriter(cm.logFile, true);
				BufferedWriter bw = new BufferedWriter(fw);
				bw.append(ds2.getID()+"\t"+k.getID()+"\t"+(explicit?"exp":"imp")+"\t"+(double)runtime/1000d/1000+" ms \n");
				bw.close();

				// write gram file 
				String fileName = cm.gramDir.getAbsolutePath()+"/"+ds2.getID()+"__"+k.getID()+".gram";
				GramUtil.writeLibSVMFile(gram, ds2.getClassLabels(), fileName);
				
				System.out.println();
			}
		}

	}
		
}
