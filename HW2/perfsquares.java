//Gaurav Anil Yeole
//EEL 6935: Distributed Computing
//Assignment 2
//UFID: 54473949

import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.Future;
import java.util.*;


class Task implements Runnable{ //All methods executed by the threads are defined in this class called Task
	private long n;
	private long k;
	
	public Task(long Num, long range){ // Constructor of class task
		this.n = Num; 
		this.k = range;
	}
	
	public long square(long N) {
		return N*N;
	}

	public long calSumSq(long n, long k){ //this method returns the sum of squares of k consecutive numbers starting from n 
		long Sum;
		Sum = sqPyramid(n+k-1) - sqPyramid(n) + square(n);
		
		return Sum;
	}

	public int isPerfectSquare(long n){ //this method returns 1 if n is perfect square
		double sqrt = Math.sqrt(n);
		int x = (int) sqrt;
		if(Math.pow(sqrt,2) == Math.pow(x,2)){
			return 1;
		}
		else 
			return 0;
	}
	
	public long sqPyramid(long n){ //this method performs the sum of squares of n consecutive numbers staring from 1
		long Sum;
		//Sum = ((n*(n+1)*((2*n)+1))/6);
		Sum = (n*(n+1));
		return ((2*n+1)*Sum)/6;
	}
	
	@Override
	public void run(){ //threads executes this method
		int perSq;
		long Sum;
		Sum = calSumSq(n,k);
		perSq = isPerfectSquare(Sum);
		if (perSq == 1){
				perfsquares.addBuffer(n);
		}
		
	}
}

public class perfsquares{
	
	public static ArrayList<Long> buffer = new ArrayList<Long>();// every thread stores number in this buffer, if sum of squares of k consecutive numbers from n is perfect square 
	
	public static synchronized void addBuffer(long num){
		buffer.add(num);
	}
	
	public static void main(String[] args){
		
		long N = Long.parseLong(args[0]);
		long k = Long.parseLong(args[1]);
		
		Runtime rt = Runtime.getRuntime();		
		int numberOfProcessors = rt.availableProcessors();
		
		ThreadPoolExecutor executor = (ThreadPoolExecutor) Executors.newFixedThreadPool(numberOfProcessors-1);
		
		for (int i=1;i<=N;i++){
			Task task = new Task(i,k);
			executor.execute(task);
		}
				
		executor.shutdown();
		while(!executor.isTerminated()){
		}
		
		if(buffer.size() == 0){
			System.out.format("There is no perfect Square such that Sum of squares of consecutive numbers for N=%d and k=%d\n",N,k);
		}
		Iterator<Long> bufferIterator = buffer.iterator();
		while(bufferIterator.hasNext()){
			System.out.println(bufferIterator.next());
		}
		
		System.out.println("All Tasks Done!!!");
	}
}