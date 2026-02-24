using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Media;
using System.Numerics;
using System.Threading;
using Ultraleap.Haptics;

namespace ContinuousSpeedExperiment
{
    public class TrialData
    {
        public int TrialNumber { get; set; }
        public int BlockNumber { get; set; }
        public float Speed { get; set; }
        public float ISI_Duration { get; set; }
        public DateTime StimulusOnset { get; set; }
    }

    public class AudioTrigger
    {
        private SoundPlayer _player;
        private bool _isEnabled;
        private string _wavFilePath;

        public AudioTrigger(string wavFilePath = "trigger_short2.wav")
        {
            _wavFilePath = wavFilePath;

            try
            {
                if (!File.Exists(_wavFilePath))
                {
                    string[] possiblePaths = new string[]
                    {
                        _wavFilePath,
                        Path.Combine(AppDomain.CurrentDomain.BaseDirectory, _wavFilePath),
                        Path.Combine(Directory.GetCurrentDirectory(), _wavFilePath),
                    };

                    bool found = false;
                    foreach (var path in possiblePaths)
                    {
                        if (File.Exists(path))
                        {
                            _wavFilePath = path;
                            found = true;
                            break;
                        }
                    }

                    if (!found)
                    {
                        Console.ForegroundColor = ConsoleColor.Yellow;
                        Console.WriteLine($"⚠ Audio trigger file not found: {wavFilePath}");
                        Console.WriteLine($"  Searched in:");
                        foreach (var path in possiblePaths)
                        {
                            Console.WriteLine($"    - {path}");
                        }
                        Console.WriteLine("  Continuing without audio triggers...");
                        Console.ResetColor();
                        _isEnabled = false;
                        return;
                    }
                }

                _player = new SoundPlayer(_wavFilePath);
                _player.LoadAsync();
                Thread.Sleep(100);

                _isEnabled = true;

                Console.ForegroundColor = ConsoleColor.Green;
                Console.WriteLine($"✓ Audio trigger loaded: {_wavFilePath}");
                Console.WriteLine($"  Full path: {Path.GetFullPath(_wavFilePath)}");
                Console.ResetColor();

                Console.WriteLine("  Testing audio... (you should hear a sound)");
                _player.Play();
                Thread.Sleep(200);
            }
            catch (Exception ex)
            {
                Console.ForegroundColor = ConsoleColor.Red;
                Console.WriteLine($"✗ Could not load audio trigger: {ex.Message}");
                Console.WriteLine($"  {ex.GetType().Name}");
                Console.ResetColor();
                _isEnabled = false;
            }
        }

        public void Play()
        {
            if (!_isEnabled) return;

            try
            {
                _player.Play();
            }
            catch (Exception ex)
            {
                Console.WriteLine($"⚠ Audio playback error: {ex.Message}");
                _isEnabled = false;
            }
        }

        public bool IsEnabled()
        {
            return _isEnabled;
        }

        public void Dispose()
        {
            _player?.Dispose();
        }
    }

    public class ContinuousSpeedExperiment
    {
        private StreamingEmitter _emitter;
        private Random _random;
        private AudioTrigger _audioTrigger;

        // Experiment parameters
        private readonly float[] speeds = { 1.0f, 3.0f, 5.0f, 7.0f, 10.0f };
        private const int TrialsPerSpeed = 150;
        private const int TrialsPerBlock = 150;
        // BLOCKS: Block 1 (1-150), Block 2 (151-300), Block 3 (301-450), Block 4 (451-600), Block 5 (601-750)

        private const float StimulusDuration = 0.2f;

        private readonly float[] possibleISIs = { 1.5f, 2.0f, 2.5f };
        private const float ISI_ChangeRate = 0.20f;

        private const float CircleDiameter = 0.065f;
        private const float CircleRadius = CircleDiameter / 2.0f;
        private readonly Vector3 _handPosition = new Vector3(0f, 0f, 0.17f);

        private List<TrialData> _allTrials;
        private StreamWriter _logFile;
        private string _csvFilename;

        public ContinuousSpeedExperiment(StreamingEmitter emitter)
        {
            _emitter = emitter;
            _random = new Random();
            _audioTrigger = new AudioTrigger("trigger_short2.wav");
        }

        public void RunExperiment(int participantId, int startBlock = 1)
        {
            Console.Clear();

            if (!_audioTrigger.IsEnabled())
            {
                Console.WriteLine("⚠ Audio triggers are NOT working!");
                Console.WriteLine("Press ENTER to continue anyway, or CTRL+C to abort...");
                Console.ReadLine();
            }

            // Generate trial sequence
            Console.WriteLine("Generating trial sequence...");
            _allTrials = GenerateTrialSequence();
            Console.WriteLine($"✓ Generated {_allTrials.Count} trials\n");

            // Write CSV before experiment starts
            _csvFilename = $"P{participantId:D3}_trials.csv";
            WriteCompleteCSV(participantId);

            InitializeLogFile(participantId, startBlock);

            int totalBlocks = (int)Math.Ceiling((double)_allTrials.Count / TrialsPerBlock);
            int startTrialIndex = (startBlock - 1) * TrialsPerBlock;
            var trialsToRun = _allTrials.Skip(startTrialIndex).ToList();

            Console.WriteLine($"Participant ID: {participantId}");
            Console.WriteLine($"Starting from: Block {startBlock}/{totalBlocks}");
            Console.WriteLine($"Trials in this session: {trialsToRun.Count} (Trial {startTrialIndex + 1} to {_allTrials.Count})");
            Console.WriteLine($"Total blocks: {totalBlocks} (150 trials each)");
            Console.WriteLine($"Duration: ~{EstimateDuration(trialsToRun.Count):F0} minutes (~{EstimateDuration(trialsToRun.Count) / 60.0f:F1} hours)");

            if (_audioTrigger.IsEnabled())
            {
                Console.WriteLine($"Audio triggers: ENABLED ✓");
            }
            else
            {
                Console.WriteLine($"Audio triggers: DISABLED ⚠");
            }

            Console.ForegroundColor = ConsoleColor.Green;
            Console.WriteLine($"✓ Complete trial list saved: {_csvFilename}");
            Console.ResetColor();
            Console.WriteLine();

            ShowInstructions();

            Console.WriteLine("\nPosition hand 17cm above device (palm down)");
            Console.WriteLine("Press ENTER when ready to start...");
            Console.ReadLine();

            LogEvent($"EXPERIMENT_START_Block{startBlock}", participantId);

            Console.Clear();
            Console.WriteLine($"Session: Block {startBlock}-{totalBlocks} | {trialsToRun.Count} trials\n");
            Console.WriteLine("Trial progress:\n");

            // Present stimuli
            foreach (var trial in trialsToRun)
            {
                // Break between blocks
                if (trial.TrialNumber > 1 &&
                    (trial.TrialNumber - 1) % TrialsPerBlock == 0 &&
                    trial.TrialNumber <= _allTrials.Count - 1)
                {
                    TakeBreak(trial.TrialNumber - 1, _allTrials.Count, trial.BlockNumber);
                }

                ExecuteTrial(trial, participantId, _allTrials.Count);
            }

            LogEvent($"EXPERIMENT_END_Block{startBlock}-{totalBlocks}", participantId);
            CloseLogFile();
            _audioTrigger.Dispose();

            Console.Clear();

            Console.WriteLine($"✓ Trial sequence saved: {_csvFilename}");
            Console.WriteLine("  (Button responses recorded by amplifier)\n");
            Console.WriteLine("\nPress any key to exit...");
            Console.ReadKey();
        }

        private void ShowInstructions()
        {
            Console.WriteLine("INSTRUCTIONS:");
            Console.WriteLine("You will feel CONTINUOUS haptic stimuli:");
            Console.WriteLine("  Stimulus → Pause → Stimulus → Pause → ...");
            Console.WriteLine();
            Console.WriteLine("YOUR TASK:");
            Console.WriteLine("  Press the button when the PAUSE duration changes");
            Console.WriteLine("  (Responses are recorded by the amplifier)");
            Console.WriteLine();
            Console.WriteLine("IMPORTANT:");
            Console.WriteLine("  • Keep your hand STILL at all times");
            Console.WriteLine($"  • Breaks after every block (150 trials)");
            Console.WriteLine($"  • This experiment has {speeds.Length * TrialsPerSpeed} trials total");
            if (_audioTrigger.IsEnabled())
            {
                Console.WriteLine("  • You will hear a sound with each stimulus");
            }
        }

        private List<TrialData> GenerateTrialSequence()
        {
            var trials = new List<TrialData>();
            int trialNum = 1;

            // Generate 150 trials for EACH speed condition
            foreach (var speed in speeds)
            {
                for (int rep = 0; rep < TrialsPerSpeed; rep++)
                {
                    trials.Add(new TrialData
                    {
                        TrialNumber = trialNum++,
                        Speed = speed
                    });
                }
            }

            // Randomize trial order
            trials = trials.OrderBy(x => _random.Next()).ToList();

            // Assign ISIs and block numbers
            float currentISI = possibleISIs[1]; // Start with middle ISI

            for (int i = 0; i < trials.Count; i++)
            {
                trials[i].TrialNumber = i + 1;
                trials[i].BlockNumber = (i / TrialsPerBlock) + 1;

                // Change ISI with 20% probability
                bool shouldChange = _random.NextDouble() < ISI_ChangeRate;

                if (shouldChange && i > 0)
                {
                    var otherISIs = possibleISIs.Where(isi => Math.Abs(isi - currentISI) > 0.1f).ToArray();
                    if (otherISIs.Length > 0)
                    {
                        currentISI = otherISIs[_random.Next(otherISIs.Length)];
                    }
                }

                trials[i].ISI_Duration = currentISI;
            }

            return trials;
        }

        // Write complete CSV
        private void WriteCompleteCSV(int participantId)
        {
            using (var writer = new StreamWriter(_csvFilename))
            {
                
                writer.WriteLine("trial,participant,block,speed,isi_duration,stimulus_onset");

                // Timestamps filled during execution
                foreach (var trial in _allTrials)
                {
                    writer.WriteLine($"{trial.TrialNumber}," +
                                   $"{participantId}," +
                                   $"{trial.BlockNumber}," +
                                   $"{trial.Speed.ToString("F1", System.Globalization.CultureInfo.InvariantCulture)}," +
                                   $"{trial.ISI_Duration.ToString("F1", System.Globalization.CultureInfo.InvariantCulture)}," +
                                   "");  // Timestamp bleibt leer
                }
            }

            Console.ForegroundColor = ConsoleColor.Green;
            Console.WriteLine($"✓ Complete trial sequence written to: {_csvFilename}");
            Console.ResetColor();
        }

        // Update CSV with timestamp during recording
        private void UpdateTrialTimestamp(TrialData trial)
        {
            try
            {
                var lines = File.ReadAllLines(_csvFilename).ToList();

                // Find line for this trial (+1 for header)
                int lineIndex = trial.TrialNumber;  ist Zeile 0, Trial 1 ist Zeile 1

                if (lineIndex < lines.Count)
                {
                    var parts = lines[lineIndex].Split(',');
                    if (parts.Length >= 6)
                    {
                        // Update timestamp column
                        parts[5] = trial.StimulusOnset.ToString("HH:mm:ss.fff");
                        lines[lineIndex] = string.Join(",", parts);

                        File.WriteAllLines(_csvFilename, lines);
                    }
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine($"⚠ Could not update CSV: {ex.Message}");
            }
        }

        private void ExecuteTrial(TrialData trial, int participantId, int totalTrials)
        {
            Console.SetCursorPosition(0, 8);

            int trialInBlock = ((trial.TrialNumber - 1) % TrialsPerBlock) + 1;
            int totalBlocks = (int)Math.Ceiling((double)totalTrials / TrialsPerBlock);

            Console.Write($"Block {trial.BlockNumber}/{totalBlocks} | Trial {trialInBlock}/{TrialsPerBlock} | ");
            Console.Write($"Overall: {trial.TrialNumber}/{totalTrials} | Speed: {trial.Speed:F1} m/s | ISI: {trial.ISI_Duration:F1}s");
            Console.Write("          ");
            Console.WriteLine();

            // Record timestamp
            trial.StimulusOnset = DateTime.Now;

            // Update CSV
            UpdateTrialTimestamp(trial);

            LogEvent($"STIM_ON_Block{trial.BlockNumber}_Speed{trial.Speed}_ISI{trial.ISI_Duration}_T{trial.TrialNumber}", participantId);

            // Audio trigger + haptic stimulus
            _audioTrigger.Play();
            PresentStimulus(trial.Speed);

            LogEvent($"STIM_OFF_T{trial.TrialNumber}", participantId);

            // ISI wait
            Thread.Sleep((int)(trial.ISI_Duration * 1000));
        }

        private void PresentStimulus(float speed)
        {
            try
            {
                DateTimeOffset startTime = DateTimeOffset.MinValue;
                bool isRunning = true;

                _emitter.EmissionCallback = (em, interval, deadline) =>
                {
                    if (!isRunning) return;
                    if (startTime == DateTimeOffset.MinValue)
                        startTime = interval.First().Time;

                    foreach (var sample in interval)
                    {
                        double elapsed = (sample.Time - startTime).TotalSeconds;
                        double phase = (speed / CircleRadius) * elapsed;

                        sample.Points[0].Position = new Vector3(
                            _handPosition.X + CircleRadius * (float)Math.Cos(phase),
                            _handPosition.Y + CircleRadius * (float)Math.Sin(phase),
                            _handPosition.Z
                        );
                        sample.Points[0].Intensity = 1.0f;
                    }
                };

                _emitter.Start();
                Thread.Sleep((int)(StimulusDuration * 1000));
                isRunning = false;
                _emitter.Stop();
                _emitter.EmissionCallback = null;
                Thread.Sleep(50);
            }
            catch (Exception ex)
            {
                Console.WriteLine($"⚠ Error: {ex.Message}");
            }
        }

        private void TakeBreak(int current, int total, int blockNumber)
        {
            int totalBlocks = (int)Math.Ceiling((double)total / TrialsPerBlock);

            Console.Clear();
            Console.WriteLine($"✓ Completed Block {blockNumber}/{totalBlocks}");
            Console.WriteLine($"✓ Completed: {current}/{total} trials ({(float)current / total * 100:F1}%)\n");

            Console.ForegroundColor = ConsoleColor.Green;
            Console.WriteLine($"✓ Trial data saved to: {_csvFilename}");
            Console.ResetColor();
            Console.WriteLine();

            Console.WriteLine($"Take a 3-5 minute break.");
            Console.WriteLine($"Remaining: {total - current} trials in {totalBlocks - blockNumber} blocks");
            Console.WriteLine();
            Console.WriteLine("Press ENTER when ready to continue...");
            Console.ReadLine();

            Console.Clear();
            Console.WriteLine($"Continuing from Block {blockNumber + 1}\n");
            Console.WriteLine("Trial progress:\n");
        }

        private void InitializeLogFile(int participantId, int startBlock)
        {
            string filename = $"P{participantId:D3}_log_Block{startBlock}_{DateTime.Now:yyyyMMdd_HHmmss}.txt";
            _logFile = new StreamWriter(filename);
            _logFile.WriteLine("Continuous Speed Experiment");
            _logFile.WriteLine($"Participant: {participantId}");
            _logFile.WriteLine($"Starting Block: {startBlock}");
            _logFile.WriteLine($"Date: {DateTime.Now}");
            _logFile.WriteLine($"Conditions: {speeds.Length} speeds");
            _logFile.WriteLine($"Trials per condition: {TrialsPerSpeed}");
            _logFile.WriteLine($"Total trials: {speeds.Length * TrialsPerSpeed}");
            _logFile.WriteLine($"Break interval: Every {TrialsPerBlock} trials");
            _logFile.WriteLine($"Audio Trigger: {(_audioTrigger.IsEnabled() ? "ENABLED" : "DISABLED")}");
            _logFile.WriteLine("Timestamp\tEvent\tParticipantID");
            _logFile.Flush();
        }

        private void LogEvent(string eventName, int participantId)
        {
            _logFile?.WriteLine($"{DateTime.Now:yyyy-MM-dd HH:mm:ss.fff}\t{eventName}\t{participantId}");
            _logFile?.Flush();
        }

        private void CloseLogFile()
        {
            _logFile?.Close();
        }

        private float EstimateDuration(int totalTrials)
        {
            float avgStimulus = StimulusDuration;
            float avgISI = possibleISIs.Average();
            float secondsPerTrial = avgStimulus + avgISI;
            int numBreaks = (totalTrials / TrialsPerBlock);
            if (totalTrials % TrialsPerBlock == 0) numBreaks--;
            float breakTime = (float)numBreaks * 180f;
            return (totalTrials * secondsPerTrial + breakTime) / 60.0f;
        }
    }

    class Program
    {
        static void Main(string[] args)
        {
            Console.OutputEncoding = System.Text.Encoding.UTF8;

            Console.Write("Connecting... ");
            Library lib = new Library();
            lib.Connect();
            IDevice device = lib.FindDevice();

            if (device == null)
            {
                Console.WriteLine("✗ No device found");
                Console.ReadKey();
                return;
            }

            StreamingEmitter emitter = new StreamingEmitter(lib);
            emitter.Devices.Add(device);
            emitter.SetControlPointCount(1, AdjustRate.None);

            // Measure device sample rate
            int totalSamples = 0;
            int totalCallbacks = 0;
            DateTimeOffset firstSampleTime = DateTimeOffset.MinValue;
            DateTimeOffset lastSampleTime = DateTimeOffset.MinValue;
            double minDt = double.MaxValue;
            double maxDt = 0;

            emitter.EmissionCallback = (em, interval, deadline) =>
            {
                foreach (var sample in interval)
                {
                    if (firstSampleTime == DateTimeOffset.MinValue)
                        firstSampleTime = sample.Time;

                    // Time between consecutive samples
                    if (lastSampleTime != DateTimeOffset.MinValue)
                    {
                        double dt = (sample.Time - lastSampleTime).TotalMilliseconds;
                        if (dt > 0)
                        {
                            if (dt < minDt) minDt = dt;
                            if (dt > maxDt) maxDt = dt;
                        }
                    }
                    lastSampleTime = sample.Time;

                    // Set dummy position (required for output)
                    sample.Points[0].Position = new Vector3(0f, 0f, 0.17f);
                    sample.Points[0].Intensity = 1.0f;

                    totalSamples++;
                }
                totalCallbacks++;
            };

            Console.WriteLine("Running 3-second sample rate test...");
            emitter.Start();
            Thread.Sleep(3000);
            emitter.Stop();
            emitter.EmissionCallback = null;

            // Evaluate results
            double totalSeconds = (lastSampleTime - firstSampleTime).TotalSeconds;
            double samplesPerSecond = totalSamples / totalSeconds;
            double samplesPerCallback = (double)totalSamples / totalCallbacks;

            Console.WriteLine($"\n═══════════════════════════════════════");
            Console.WriteLine($"  RESULTS:");
            Console.WriteLine($"  Total duration:          {totalSeconds:F3} s");
            Console.WriteLine($"  Total samples:       {totalSamples}");
            Console.WriteLine($"  Total callbacks:     {totalCallbacks}");
            Console.WriteLine($"  ─────────────────────────────────────");
            Console.WriteLine($"  Sample rate:          {samplesPerSecond:F0} Hz");
            Console.WriteLine($"  Samples per callback: {samplesPerCallback:F1}");
            Console.WriteLine($"  Min dt (Samples):     {minDt:F4} ms");
            Console.WriteLine($"  Max dt (Samples):     {maxDt:F4} ms");

            Thread.Sleep(100);
            Console.WriteLine("✓\n");

            Console.Write("Enter Participant ID: ");
            int participantId;
            if (!int.TryParse(Console.ReadLine(), out participantId))
            {
                participantId = 999;
            }

            Console.Write("\nStart from Block (1-5, default=1): ");
            string blockInput = Console.ReadLine();
            int startBlock = 1;
            if (!string.IsNullOrWhiteSpace(blockInput) && int.TryParse(blockInput, out int inputBlock))
            {
                if (inputBlock >= 1 && inputBlock <= 5)
                {
                    startBlock = inputBlock;
                }
                else
                {
                    Console.WriteLine("Invalid block number. Starting from Block 1.");
                }
            }

            try
            {
                var experiment = new ContinuousSpeedExperiment(emitter);
                experiment.RunExperiment(participantId, startBlock);
            }
            catch (Exception ex)
            {
                Console.WriteLine($"\n❌ ERROR: {ex.Message}");
                Console.WriteLine($"\n{ex.StackTrace}");
            }
            finally
            {
                emitter.Stop();
                emitter.Devices.Remove(device);
            }

            Console.ReadKey();
        }
    }
}