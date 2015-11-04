using System;
using System.IO;
using System.Threading;
using System.Text;
using System.Net;
using System.Net.Sockets;
using System.Collections;
using System.Text.RegularExpressions;

namespace My
{
    public class NetThread
    {
		// Receiving Thread Object
		Thread m_receiveThread;
		
		// UDP Client Object
		UDPClient m_client;

		// IPEndPoint
		IPEndPoint endPoint;
		
		// IP Address and Port
		private string IP = "127.0.0.1";	// default local
		private int port; 					// define in init
		
		// Containers for the recieved UDP packets
        private float[,] oldPoints = new float[26, 3] { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, 
                 { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
                 { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };

        private float[,] newPoints = new float[26, 3] { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, 
                 { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
                 { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };

        // Container for the marker short codes             Thumb                       Index                    Middle                 Ring                  Little                
        private string[] shortCodes = new string[26] { "T0", "T1", "T2", "T3", "I0", "I1", "I2", "I3", "M0", "M1", "M2", "M3", "R0", "R1", "R2", "R3", "L0", "L1", "L2", "L3",
                                                       "HB", "WJ", "FA",    // Hand Back Center, Wrist joint, Fore Arm,  
                                                       "AD", "AL", "AV"};   // Amputation Distal, Amputation lateral, Amputation Ventral/Dorsal

		// Volatile bool for termination of thread
		public volatile bool terminateFlag;
		
		// String Delimeters
		char[] delimiterChars;
			
	    // Creates the pool of semaphores
	    private static Semaphore pool;

		// Start from Matlab
		void Start(){
			init();
		}
		
		// Initialization of variables
		public void init(){
    		// define port
			port = 9080;
			
			//		delimiterChars = new char[] {':', ';', '\t', '\n'};
			delimiterChars = new char[] {':', ';'};
			
			// Set termination boolean
			terminateFlag = false;

		    // Creates a new semaphore pool with only 1 thread availiable
		    pool = new Semaphore(0, 1);
		    pool.Release(1);
			
		}// END FUNCTION

        // Set IP address
        public void SetIP(string tmpIP) {
            IP = tmpIP;
        } // END FUNCTION SetIP

        // Set IP port
        public void SetPort(int tmpPort) {
            port = tmpPort;
        } // END FUNCTION SetPort

	    public void CreateThread(){
	    	// Create thread
	        m_receiveThread = new Thread(new ThreadStart(ReceiveData)); // Set function to thread
	    	m_receiveThread.IsBackground = true;						// Run in background
	    	m_receiveThread.Start();	
	    }// END FUNCTION

        // Set termination flag to foce the thread to shutdown
        public void SetTerminateFlag(bool tmpBool)
        { terminateFlag = tmpBool; } // END FUNCTION

        // Set the latest points from the Vicon system
        private void SetPoints(float[,] tmpPoints) {
            pool.WaitOne();
            oldPoints = tmpPoints;
            pool.Release();
        } // END FUNCTION SetPoints()

        // Get the letest points from the Vicon system
        public float[,] GetPoints() {
            pool.WaitOne();
            float[,] tmpPoints = oldPoints;
            pool.Release();
            return tmpPoints;
        } // END FUNCTION GetPoints()

        // Set the short codes used in the Vicon packet streaming
        public void SetCodes(string[] tmpString) {
            shortCodes = tmpString;
        } // END FUNCTION SetCodes()

        // Parse the Vicon UDP packet on a separate thread
        private void ReciveveData() {
            try {
                while(!this.terminateFlag) {

                    IPEndPoint endPoint = new IPEndPoint(IPAddress.Parse(IP), port);
                    byte[] data = m_client.Recieve(ref endPoint);
                        
                    // Convert bytes to UTF-8 string
                    string text = Encoding.UTF8.GetString(data);

                    // Initialize a verbatium expression string
                    string expr = @"";

                    // Initialize a container for float parsing
                    float tmpNum = 0;

                    // For each short code, regex through the UDP packet
                    for (int ii = 0; ii < shortCodes.Length; ii++){

                        // Grab the short code
                        string code = shortCodes[ii];

                        // Build expression string
                        expr = code + @":\t([0-9\.]+)\t([0-9\.]+)\t([0-9\.]+)";

                        // Attempt regex
                        regOut = regex.Match(text, expr);

                        // IF the regex found a match, then parse x-y-z
                        if (regOut != NULL) {
                            // Attempt to parse the x-value as a float
                            bool resX = float.TryParse(regOut[1], out tempNum);
                            if (resX == true) {
                                newPoints[ii, 1] = tempNum;
                            } // END IF TryParse(x)

                            // Attempt to parse the y-value as a float
                            bool resY = float.TryParse(regOut[2], out tempNum);
                            if (resY == true){
                                newPoints[ii, 2] = tempNum;
                            } // END IF TryParse(y)

                            // Attempt to parse the z-value as a float
                            bool resZ = float.TryParse(regOut[3], out tempNum);
                            if (resZ == true) {
                                newPoints[ii, 3] = tempNum;
                            } // END IF TryParse(z)
                        } // END IF regex output
                    } // END FOR shortCodes

                    // Set newPoints
                    this.SetPoints(newPoints);

                } // END FOR
            } // END TRY
            catch (Exception ex) {
                Console.WriteLine(ex.Message); // Output error message
            } // END CATCH
        } // END FUNCTION
    } // END CLASS
} // END NAMESPACE

// EOF