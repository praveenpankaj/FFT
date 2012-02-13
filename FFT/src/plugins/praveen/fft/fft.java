package plugins.praveen.fft;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_2D;
import edu.emory.mathcs.jtransforms.fft.DoubleFFT_3D;
import icy.image.IcyBufferedImage;
import icy.sequence.Sequence;
import icy.type.DataType;
import icy.type.collection.array.Array1DUtil;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.ezplug.EzVarText;

public class fft extends EzPlug {

	EzVarSequence input = new EzVarSequence("Input");
	EzVarText	ndims = new EzVarText("Type", new String[] { "2D", "3D" }, 0, false);
	EzVarText	display = new EzVarText("Display as", new String[] {  "Magnitude/Phase Pair", "Real/Imaginary Pair" }, 0, false);
	EzVarText	swap = new EzVarText("Swap Quadrants?", new String[] { "Yes", "No" }, 1, false);

	@Override
	protected void initialize() {
		super.addEzComponent(input);
		super.addEzComponent(ndims);
		super.addEzComponent(swap);
		super.addEzComponent(display);		
		super.setTimeDisplay(true);
	}

	@Override
	protected void execute() {
		Sequence sequence = input.getValue();

		if(ndims.getValue()=="2D")		
			FFT_2D(sequence, swap.getValue(), display.getValue());	
		else
			FFT_3D(sequence, swap.getValue(), display.getValue());
		//MessageDialog.showDialog("FFT3D not implemented yet !");	
	}



	private Sequence FFT_3D(Sequence sequence, String swap, String display) {
		int _w = sequence.getSizeX();
		int _h = sequence.getSizeY();
		int _z = sequence.getSizeZ();
		int wc = (int) Math.ceil(_w/2);
		int hc = (int) Math.ceil(_h/2);
		int zc = (int) Math.ceil(_z/2);

		double[] fArray;
		final DoubleFFT_3D fft = new DoubleFFT_3D(_z, _h, _w);
		Sequence fSequence = new Sequence();
		fSequence.setName("Fourier Transform 3D");

		if(swap == "No")
		{ //No Quadrant swapping. Leave as it is.
			for(int k = 0; k < _z; k++)
			{	
				IcyBufferedImage resultMatrix = new IcyBufferedImage(_w, _h, 2, DataType.DOUBLE);			
				resultMatrix.setDataXY(0, Array1DUtil.arrayToDoubleArray(sequence.getDataXY(0, k, 0), sequence.isSignedDataType()));//set buffered image to sequence 
				fSequence.setImage(0, k, resultMatrix);			
			}						
			fArray = fSequence.getDataCopyCXYZAsDouble(0);
			fft.complexForward(fArray);//Does only on half the data. To get the full transform use realForwardFull

			if(display=="Magnitude/Phase Pair")
			{
				for(int k = 0; k < _z; k++)
				{			
					IcyBufferedImage resultMatrix = fSequence.getImage(0, k);


					resultMatrix.beginUpdate();
					try
					{
						for(int x = 0; x < _w; x++)
						{
							for(int y = 0; y < _h; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(fArray[(x + (y * _w) + (k * _w * _h))*2 + 0], 2)+Math.pow(fArray[(x + (y * _w) + (k * _w * _h))*2 + 1], 2)));
								resultMatrix.setDataAsDouble(x, y, 1, Math.atan2(fArray[(x + (y * _w) + (k * _w * _h))*2 + 1], fArray[(x + (y * _w) + (k * _w * _h))*2 + 0]));
							}
						}
					}
					finally
					{
						resultMatrix.endUpdate();
					}
				}
			}
			else
			{
				for(int k = 0; k < _z; k++)
				{			
					IcyBufferedImage resultMatrix = fSequence.getImage(0, k);
					resultMatrix.beginUpdate();
					try
					{
						for(int x = 0; x < _w; x++)
						{
							for(int y = 0; y < _h; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, fArray[(x + (y * _w) + (k * _w * _h))*2 + 0]);
								resultMatrix.setDataAsDouble(x, y, 1, fArray[(x + (y * _w) + (k * _w * _h))*2 + 1]);
							}
						}
					}
					finally
					{
						resultMatrix.endUpdate();
					}
				}
			}
		}
		else
		{//Swap Quadrants
			for(int k = 0; k < _z; k++)
			{	
				IcyBufferedImage resultMatrix = new IcyBufferedImage(_w, _h, 2, DataType.DOUBLE);			
				resultMatrix.setDataXY(0, Array1DUtil.arrayToDoubleArray(sequence.getDataXY(0, k, 0), sequence.isSignedDataType()));//set buffered image to sequence 
				fSequence.setImage(0, k, resultMatrix);			
			}						
			fArray = fSequence.getDataCopyCXYZAsDouble(0);
			fft.complexForward(fArray);//Does only on half the data. To get the full transform use realForwardFull

			if(display=="Magnitude/Phase Pair")
			{
				for(int k = 0; k < zc+1; k++)
				{			
					IcyBufferedImage resultMatrix = fSequence.getImage(0, k);


					resultMatrix.beginUpdate();
					try
					{
						for(int x = 0; x < wc+1; x++)
						{
							for(int y = 0; y < hc+1; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(fArray[((wc-x) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 0], 2)+Math.pow(fArray[((wc-x) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 1], 2)));
								resultMatrix.setDataAsDouble(x, y, 1, Math.atan2(fArray[((wc-x) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 1], fArray[((wc-x) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 0]));
							}
							for(int y = hc+1; y < _h; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(fArray[((wc-x) + (y-hc) * _w + (zc-k) * _w * _h)*2 + 0], 2)+Math.pow(fArray[((wc-x) + (y-hc) * _w + (zc-k) * _w * _h)*2 + 1], 2)));
								resultMatrix.setDataAsDouble(x, y, 1, Math.atan2(fArray[((wc-x) + (y-hc) * _w + (zc-k) * _w * _h)*2 + 1], fArray[((wc-x) + (y-hc) * _w + (zc-k) * _w * _h)*2 + 0]));
							}
						}
						for(int x = wc+1; x < _w; x++)
						{
							for(int y = 0; y < hc+1; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(fArray[((x-wc) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 0], 2)+Math.pow(fArray[((x-wc) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 1], 2)));
								resultMatrix.setDataAsDouble(x, y, 1, Math.atan2(fArray[((x-wc) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 1], fArray[((x-wc) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 0]));
							}
							for(int y = hc+1; y < _h; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(fArray[((x-wc) + (y-hc) * _w + (zc-k) * _w * _h)*2 + 0], 2)+Math.pow(fArray[((x-wc) + (y-hc) * _w + (zc-k) * _w * _h)*2 + 1], 2)));
								resultMatrix.setDataAsDouble(x, y, 1, Math.atan2(fArray[((x-wc) + (y-hc) * _w + (zc-k) * _w * _h)*2 + 1], fArray[((x-wc) + (y-hc) * _w + (zc-k) * _w * _h)*2 + 0]));
							}
						}
					}

					finally
					{
						resultMatrix.endUpdate();
					}
				}
				for(int k = zc+1; k < _z; k++)
				{			
					IcyBufferedImage resultMatrix = fSequence.getImage(0, k);


					resultMatrix.beginUpdate();
					try
					{
						for(int x = 0; x < wc+1; x++)
						{
							for(int y = 0; y < hc+1; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(fArray[((wc-x) + (hc-y) * _w + (k-zc) * _w * _h)*2 + 0], 2)+Math.pow(fArray[((wc-x) + (hc-y) * _w + (k-zc) * _w * _h)*2 + 1], 2)));
								resultMatrix.setDataAsDouble(x, y, 1, Math.atan2(fArray[((wc-x) + (hc-y) * _w + (k-zc) * _w * _h)*2 + 1], fArray[((wc-x) + (hc-y) * _w + (k-zc) * _w * _h)*2 + 0]));
							}
							for(int y = hc+1; y < _h; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(fArray[((wc-x) + (y-hc) * _w + (k-zc) * _w * _h)*2 + 0], 2)+Math.pow(fArray[((wc-x) + (y-hc) * _w + (k-zc) * _w * _h)*2 + 1], 2)));
								resultMatrix.setDataAsDouble(x, y, 1, Math.atan2(fArray[((wc-x) + (y-hc) * _w + (k-zc) * _w * _h)*2 + 1], fArray[((wc-x) + (y-hc) * _w + (k-zc) * _w * _h)*2 + 0]));
							}
						}
						for(int x = wc+1; x < _w; x++)
						{
							for(int y = 0; y < hc+1; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(fArray[((x-wc) + (hc-y) * _w + (k-zc) * _w * _h)*2 + 0], 2)+Math.pow(fArray[((x-wc) + (hc-y) * _w + (k-zc) * _w * _h)*2 + 1], 2)));
								resultMatrix.setDataAsDouble(x, y, 1, Math.atan2(fArray[((x-wc) + (hc-y) * _w + (k-zc) * _w * _h)*2 + 1], fArray[((x-wc) + (hc-y) * _w + (k-zc) * _w * _h)*2 + 0]));
							}
							for(int y = hc+1; y < _h; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(fArray[((x-wc) + (y-hc) * _w + (k-zc) * _w * _h)*2 + 0], 2)+Math.pow(fArray[((x-wc) + (y-hc) * _w + (k-zc) * _w * _h)*2 + 1], 2)));
								resultMatrix.setDataAsDouble(x, y, 1, Math.atan2(fArray[((x-wc) + (y-hc) * _w + (k-zc) * _w * _h)*2 + 1], fArray[((x-wc) + (y-hc) * _w + (k-zc) * _w * _h)*2 + 0]));
							}
						}
					}
					finally
					{
						resultMatrix.endUpdate();
					}
				}
			}
			else
			{
				for(int k = 0; k < zc+1; k++)
				{			
					IcyBufferedImage resultMatrix = fSequence.getImage(0, k);
					resultMatrix.beginUpdate();
					try
					{
						for(int x = 0; x < wc+1; x++)
						{
							for(int y = 0; y < hc+1; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, fArray[((wc-x) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 0]);
								resultMatrix.setDataAsDouble(x, y, 1, fArray[((wc-x) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 1]);
							}
							for(int y = hc+1; y < _h; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, fArray[((wc-x) + (y-hc) * _w + (zc-k) * _w * _h)*2 + 0]);
								resultMatrix.setDataAsDouble(x, y, 1, fArray[((wc-x) + (y-hc) * _w + (zc-k) * _w * _h)*2 + 1]);
							}
						}
						for(int x = wc+1; x < _w; x++)
						{
							for(int y = 0; y < hc+1; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, fArray[((x-wc) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 0]);
								resultMatrix.setDataAsDouble(x, y, 1, fArray[((x-wc) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 1]);
							}
							for(int y = hc+1; y < _h; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, fArray[((x-wc) + (y-hc) * _w + (zc-k) * _w * _h)*2 + 0]);
								resultMatrix.setDataAsDouble(x, y, 1, fArray[((x-wc) + (y-hc) * _w + (zc-k) * _w * _h)*2 + 1]);
							}
						}
					}
					finally
					{
						resultMatrix.endUpdate();
					}
				}
				for(int k = zc+1; k < _z; k++)
				{			
					IcyBufferedImage resultMatrix = fSequence.getImage(0, k);
					resultMatrix.beginUpdate();
					try
					{
						for(int x = 0; x < wc+1; x++)
						{
							for(int y = 0; y < hc+1; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, fArray[((wc-x) + (hc-y) * _w + (k-zc) * _w * _h)*2 + 0]);
								resultMatrix.setDataAsDouble(x, y, 1, fArray[((wc-x) + (hc-y) * _w + (k-zc) * _w * _h)*2 + 1]);
							}
							for(int y = hc+1; y < _h; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, fArray[((wc-x) + (y-hc) * _w + (k-zc) * _w * _h)*2 + 0]);
								resultMatrix.setDataAsDouble(x, y, 1, fArray[((wc-x) + (y-hc) * _w + (k-zc) * _w * _h)*2 + 1]);
							}
						}
						for(int x = wc+1; x < _w; x++)
						{
							for(int y = 0; y < hc+1; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, fArray[((x-wc) + (hc-y) * _w + (k-zc) * _w * _h)*2 + 0]);
								resultMatrix.setDataAsDouble(x, y, 1, fArray[((x-wc) + (hc-y) * _w + (k-zc) * _w * _h)*2 + 1]);
							}
							for(int y = hc+1; y < _h; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, fArray[((x-wc) + (y-hc) * _w + (k-zc) * _w * _h)*2 + 0]);
								resultMatrix.setDataAsDouble(x, y, 1, fArray[((x-wc) + (y-hc) * _w + (k-zc) * _w * _h)*2 + 1]);
							}
						}
					}
					finally
					{
						resultMatrix.endUpdate();
					}
				}

			}

		}
		addSequence(fSequence);
		return fSequence;
	}

	private Sequence FFT_2D(Sequence sequence, String swap, String display) 
	{
		Sequence fSequence = new Sequence();
		fSequence.setName("Fourier Transform 2D");
		int _w = sequence.getSizeX();
		int _h = sequence.getSizeY();
		int _z = sequence.getSizeZ();
		int wc = (int) Math.ceil(_w/2);
		int hc = (int) Math.ceil(_h/2);

		final DoubleFFT_2D fft = new DoubleFFT_2D(_h, _w);
		if(swap == "No")
		{ //No Quadrant swapping

			for(int k = 0; k < _z; k++)
			{
				IcyBufferedImage resultArray = new IcyBufferedImage(_w, _h, 2, DataType.DOUBLE);

				resultArray.setDataXY(0, Array1DUtil.arrayToDoubleArray(sequence.getDataXY(0, k, 0), sequence.isSignedDataType()));
				double[] fArray = resultArray.getDataCopyCXYAsDouble();
				fft.complexForward(fArray);//Does only on half the data. To get the full transform use realForwardFull
				resultArray.beginUpdate();
				try
				{
					if(display=="Magnitude/Phase Pair")
					{
						fSequence.setChannelName(0, "Magnitude");
						fSequence.setChannelName(1, "Phase");
						for(int x = 0; x < _w; x++)
						{
							for(int y = 0; y < _h; y++)
							{
								resultArray.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(fArray[(x + y * _w)*2 + 0], 2) + Math.pow(fArray[(x + y * _w)*2 + 1], 2)));
								resultArray.setDataAsDouble(x, y, 1, Math.atan2(fArray[(x + y * _w)*2 + 1], fArray[(x + y * _w)*2 + 0]));
							}
						}

						//fImage.setChannelName(0, "Magnitude");
						//fImage.setChannelName(1, "Phase");
					}
					else // Real/Imaginary Pair
					{	fSequence.setChannelName(0, "Real");
					fSequence.setChannelName(1, "Imaginary");
						for(int x = 0; x < _w; x++)
						{
							for(int y = 0; y < _h; y++)
							{
								resultArray.setDataAsDouble(x, y, 0, fArray[(x + y * _w)*2 + 0]);
								resultArray.setDataAsDouble(x, y, 1, fArray[(x + y * _w)*2 + 1]);
							}
						}

					}

				}finally{
					resultArray.endUpdate();
				}

				fSequence.setImage(0, k, resultArray);
			}
			
		}
		else
		{ //Swap quadrants
			for(int k = 0; k < _z; k++)
			{
				IcyBufferedImage resultArray = new IcyBufferedImage(_w, _h, 2, DataType.DOUBLE);

				resultArray.setDataXY(0, Array1DUtil.arrayToDoubleArray(sequence.getDataXY(0, k, 0), sequence.isSignedDataType()));
				double[] fArray = resultArray.getDataCopyCXYAsDouble();
				fft.complexForward(fArray);//Does only on half the data. To get the full transform use realForwardFull
				resultArray.beginUpdate();
				try
				{
					if(display=="Magnitude/Phase Pair")
					{fSequence.setChannelName(0, "Magnitude");
					fSequence.setChannelName(1, "Phase");
						for(int x = 0; x < (wc+1); x++)
						{
							for(int y = 0; y < (hc+1); y++)
							{
								resultArray.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(fArray[((wc-x) + (hc-y) * _w)*2 + 0], 2)+Math.pow(fArray[((wc-x) + (hc-y) * _w)*2 + 1], 2)));
								resultArray.setDataAsDouble(x, y, 1, Math.atan2(fArray[((wc-x) + (hc-y) * _w)*2 + 1], fArray[((wc-x) + (hc-y) * _w)*2 + 0]));
							}
							for(int y = hc+1; y < _h; y++)
							{
								resultArray.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(fArray[((wc-x) + (y-hc) * _w)*2 + 0], 2)+Math.pow(fArray[((wc-x) + (y-hc) * _w)*2 + 1], 2)));
								resultArray.setDataAsDouble(x, y, 1, Math.atan2(fArray[((wc-x) + (y-hc) * _w)*2 + 1], fArray[((wc-x) + (y-hc) * _w)*2 + 0]));
							}

						}
						for(int x = (wc+1); x < _w; x++)
						{
							for(int y = 0; y < (hc+1); y++)
							{
								resultArray.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(fArray[((x-wc) + (hc-y) * _w)*2 + 0], 2)+Math.pow(fArray[((x-wc) + (hc-y) * _w)*2 + 1], 2)));
								resultArray.setDataAsDouble(x, y, 1, Math.atan2(fArray[((x-wc) + (hc-y) * _w)*2 + 1], fArray[((x-wc) + (hc-y) * _w)*2 + 0]));
							}
							for(int y = hc+1; y < _h; y++)
							{
								resultArray.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(fArray[((x-wc) + (y-hc) * _w)*2 + 0], 2)+Math.pow(fArray[((x-wc) + (y-hc) * _w)*2 + 1], 2)));
								resultArray.setDataAsDouble(x, y, 1, Math.atan2(fArray[((x-wc) + (y-hc) * _w)*2 + 1], fArray[((x-wc) + (y-hc) * _w)*2 + 0]));
							}
						}
					}					
					else //Real/Imaginary Pair
					{fSequence.setChannelName(0, "Real");
					fSequence.setChannelName(1, "Imaginary");
						for(int x = 0; x < (wc+1); x++)
						{
							for(int y = 0; y < (hc+1); y++)
							{
								resultArray.setDataAsDouble(x, y, 0, fArray[((wc-x) + (hc-y) * _w)*2 + 0]);
								resultArray.setDataAsDouble(x, y, 1, fArray[((wc-x) + (hc-y) * _w)*2 + 1]);
							}
							for(int y = hc+1; y < _h; y++)
							{
								resultArray.setDataAsDouble(x, y, 0, fArray[((wc-x) + (y-hc) * _w)*2 + 0]);
								resultArray.setDataAsDouble(x, y, 1, fArray[((wc-x) + (y-hc) * _w)*2 + 1]);
							}

						}
						for(int x = (wc+1); x < _w; x++)
						{
							for(int y = 0; y < (hc+1); y++)
							{
								resultArray.setDataAsDouble(x, y, 0, fArray[((x-wc) + (hc-y) * _w)*2 + 0]);
								resultArray.setDataAsDouble(x, y, 1, fArray[((x-wc) + (hc-y) * _w)*2 + 1]);
							}
							for(int y = hc+1; y < _h; y++)
							{
								resultArray.setDataAsDouble(x, y, 0, fArray[((x-wc) + (y-hc) * _w)*2 + 0]);
								resultArray.setDataAsDouble(x, y, 1, fArray[((x-wc) + (y-hc) * _w)*2 + 1]);
							}
						}						

					}

				}finally{
					resultArray.endUpdate();
				}

				fSequence.setImage(0, k, resultArray);
			}
			
		}


		addSequence(fSequence);

		return fSequence;
	}

	@Override
	public void clean() {
	}
}