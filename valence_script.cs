using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Threading;
using TcpServer_And_Ultrahaptics;
using Ultraleap.Haptics;


namespace ContinuousSpeedExperiment
{
    public class Trial
    {
        public int Number, Block;
        public float SpeedFirst, SpeedSecond;
        public int Choice;
        public float ChosenSpeed;
        public double RT;
        public DateTime Timestamp;
    }

    public class PairedComparisonExperiment
    {
        private StreamingEmitter _emitter;
        private Random _rng = new Random();

        private readonly float[] _speeds = { 1.0f, 3.0f, 5.0f, 7.0f, 10.0f };
        private const int RepsPerDirection = 10;
        private const int BlockSize = 40;

        private const float StimDuration = 0.2f;
        private const float CircleRadius = 0.065f / 2.0f;
        private const int ISI_ms = 1000;
        private const int ITI_ms = 800;
        private readonly Vector3 _handPos = new Vector3(0f, 0f, 0.17f);

        private List<Trial> _trials;
        private List<Trial> _done = new List<Trial>();
        private string _csv;
        private DateTime _startTime;

        public PairedComparisonExperiment(StreamingEmitter emitter)
        {
            _emitter = emitter;
        }

        public void Run(int participantId)
        {
            _trials = GenerateTrials();
            int totalBlocks = (int)Math.Ceiling((double)_trials.Count / BlockSize);

            Console.Clear();
            PrintHeader("PAIRED COMPARISON — VALENCE");
            Console.WriteLine($"  {_trials.Count} Trials ({_speeds.Length * (_speeds.Length - 1) / 2} pairs x 2 directions x {RepsPerDirection} reps)");
            Console.WriteLine($"  {totalBlocks} blocks of {BlockSize} trials, break between blocks");
            Console.WriteLine($"  Question: \"Which felt MORE PLEASANT?\"");
            Console.WriteLine($"  Stimuli werden automatisch praesentiert\n");

            _csv = $"P{participantId:D3}_valence_{DateTime.Now:yyyyMMdd_HHmmss}.csv";
            using (var w = new StreamWriter(_csv, false))
                w.WriteLine("trial,block,participant,speed_first,speed_second,choice,chosen_speed,rt_ms,timestamp");

            Console.WriteLine($"  CSV: {_csv}\n");
            Console.WriteLine("  CONTROLS:");
            Console.WriteLine("    1   = Choose Stimulus 1");
            Console.WriteLine("    0   = No clear preference");
            Console.WriteLine("    2   = Choose Stimulus 2");
            Console.WriteLine("    Esc = Abbrechen\n");
            Console.WriteLine("  Position hand 17 cm above the device.");
            Console.WriteLine("  Press ENTER to start...");
            Console.ReadLine();

            _startTime = DateTime.Now;

            for (int i = 0; i < _trials.Count; i++)
            {
                if (i > 0 && i % BlockSize == 0)
                    TakeBreak(i, _trials.Count);

                var t = _trials[i];
                ShowStatus(t);
                PresentPair(t);

                Console.ForegroundColor = ConsoleColor.Green;
                Console.WriteLine("  >>> Which felt MORE PLEASANT?  [1] or [2]");
                Console.ResetColor();

                var (ch, rt) = WaitForChoice();
                if (ch == 0) { Finish(true); return; }

                t.Choice = ch;
                t.ChosenSpeed = ch == -1 ? 0f : (ch == 1 ? t.SpeedFirst : t.SpeedSecond);
                t.RT = rt;
                t.Timestamp = DateTime.Now;

                Console.ForegroundColor = ConsoleColor.DarkGray;
                Console.WriteLine(ch == -1 ? $"    -> No preference ({rt:F0} ms)" : $"    -> Stim {ch} ({t.ChosenSpeed:F1} m/s) {rt:F0} ms");
                Console.ResetColor();

                _done.Add(t);
                AppendRow(t, participantId);

                Thread.Sleep(ITI_ms);
            }

            Finish(false);
        }

        private void ShowStatus(Trial t)
        {
            Console.Clear();
            int totalBlocks = (int)Math.Ceiling((double)_trials.Count / BlockSize);
            double elapsed = (DateTime.Now - _startTime).TotalMinutes;

            Console.ForegroundColor = ConsoleColor.Cyan;
            Console.WriteLine($"  Trial {t.Number}/{_trials.Count}  |  Block {t.Block}/{totalBlocks}  |  {elapsed:F0} min");
            Console.ResetColor();
            Console.ForegroundColor = ConsoleColor.DarkGray;
            Console.WriteLine($"  [Stim1 = {t.SpeedFirst:F1} m/s  |  Stim2 = {t.SpeedSecond:F1} m/s]");
            Console.ResetColor();
            Console.ForegroundColor = ConsoleColor.Green;
            Console.WriteLine($"\n  VALENCE — Which is more pleasant?\n");
            Console.ResetColor();
        }

        private void PresentPair(Trial t)
        {
            Console.Write("  >>> Stimulus 1 ...");
            PlayStim(t.SpeedFirst);
            Console.WriteLine(" ok");

            Console.Write("  --- Pause ---");
            Thread.Sleep(ISI_ms);
            Console.WriteLine();

            Console.Write("  >>> Stimulus 2 ...");
            PlayStim(t.SpeedSecond);
            Console.WriteLine(" ok\n");
        }

        private void PlayStim(float speed)
        {
            DateTimeOffset t0 = DateTimeOffset.MinValue;
            bool go = true;

            _emitter.EmissionCallback = (em, iv, dl) =>
            {
                if (!go) return;
                if (t0 == DateTimeOffset.MinValue) t0 = iv.First().Time;
                foreach (var s in iv)
                {
                    double el = (s.Time - t0).TotalSeconds;
                    double ph = (speed / CircleRadius) * el;
                    s.Points[0].Position = new Vector3(
                        _handPos.X + CircleRadius * (float)Math.Cos(ph),
                        _handPos.Y + CircleRadius * (float)Math.Sin(ph),
                        _handPos.Z);
                    s.Points[0].Intensity = 1.0f;
                }
            };

            _emitter.Start();
            Thread.Sleep((int)(StimDuration * 1000));
            go = false;
            _emitter.Stop();
            _emitter.EmissionCallback = null;
            Thread.Sleep(50);
        }

        private (int, double) WaitForChoice()
        {
            var t = DateTime.Now;
            while (true)
            {
                var k = Console.ReadKey(true);
                double rt = (DateTime.Now - t).TotalMilliseconds;
                if (k.Key == ConsoleKey.D0 || k.Key == ConsoleKey.NumPad0) return (-1, rt);
                if (k.Key == ConsoleKey.D1 || k.Key == ConsoleKey.NumPad1) return (1, rt);
                if (k.Key == ConsoleKey.D2 || k.Key == ConsoleKey.NumPad2) return (2, rt);
                if (k.Key == ConsoleKey.Escape) return (0, 0);
            }
        }

        private void TakeBreak(int done, int total)
        {
            double min = (DateTime.Now - _startTime).TotalMinutes;
            Console.Clear();
            Console.WriteLine($"  Break — {done}/{total} trials ({100.0 * done / total:F0}%) | {min:F0} min elapsed");
            Console.WriteLine("  Press ENTER to continue...");
            Console.ReadLine();
        }

        private void Finish(bool aborted)
        {
            double min = (DateTime.Now - _startTime).TotalMinutes;
            Console.Clear();
            PrintHeader(aborted ? "ABORTED" : "EXPERIMENT COMPLETE");
            Console.WriteLine($"  {_done.Count} Trials in {min:F1} min");
            Console.WriteLine($"  CSV: {_csv}\n");

            if (_done.Count > 0)
            {
                Console.WriteLine("  Win rate per speed:");
                foreach (var sp in _speeds)
                {
                    int w = _done.Count(t => Math.Abs(t.ChosenSpeed - sp) < 0.01f);
                    int a = _done.Count(t => Math.Abs(t.SpeedFirst - sp) < 0.01f || Math.Abs(t.SpeedSecond - sp) < 0.01f);
                    Console.WriteLine($"    {sp,5:F1} m/s: {w,3}/{a,3} ({(a > 0 ? 100.0 * w / a : 0):F1}%)");
                }
                int n = _done.Count;
                int noOp = _done.Count(t => t.Choice == -1);
                int c1 = _done.Count(t => t.Choice == 1);
                Console.WriteLine($"\n  Position: Stim1={c1} Stim2={n - c1 - noOp} No preference={noOp} ({100.0 * noOp / n:F1}%)");
                Console.WriteLine($"  RT: {_done.Average(t => t.RT):F0}ms");
            }

            Console.WriteLine("\n  Press any key to exit...");
            Console.ReadKey();
        }

        private List<Trial> GenerateTrials()
        {
            var list = new List<Trial>();
            for (int r = 0; r < RepsPerDirection; r++)
                for (int i = 0; i < _speeds.Length; i++)
                    for (int j = i + 1; j < _speeds.Length; j++)
                    {
                        list.Add(new Trial { SpeedFirst = _speeds[i], SpeedSecond = _speeds[j] });
                        list.Add(new Trial { SpeedFirst = _speeds[j], SpeedSecond = _speeds[i] });
                    }
            for (int k = list.Count - 1; k > 0; k--)
            { int s = _rng.Next(k + 1); var tmp = list[k]; list[k] = list[s]; list[s] = tmp; }
            for (int k = 0; k < list.Count; k++)
            { list[k].Number = k + 1; list[k].Block = k / BlockSize + 1; }
            return list;
        }

        private void AppendRow(Trial t, int pid)
        {
            var ci = CultureInfo.InvariantCulture;
            using (var w = new StreamWriter(_csv, true))
                w.WriteLine(string.Join(",", t.Number, t.Block, $"P{pid:D3}",
                    t.SpeedFirst.ToString("F1", ci), t.SpeedSecond.ToString("F1", ci),
                    t.Choice, t.ChosenSpeed.ToString("F1", ci), t.RT.ToString("F0", ci),
                    t.Timestamp.ToString("yyyy-MM-dd HH:mm:ss.fff")));
        }

        static void PrintHeader(string title) => Console.WriteLine($"\n{title}\n{new string('-', title.Length)}\n");
    }

    class Program
    {
        static void Main(string[] args)
        {
            Console.OutputEncoding = System.Text.Encoding.UTF8;
            PrintHeader("PAIRED COMPARISON — VALENCE");

            Console.Write("  Connecting... ");
            Library lib = new Library();
            lib.Connect();
            IDevice device = lib.FindDevice();

            if (device == null) { Console.WriteLine("No device found."); Console.ReadKey(); return; }

            StreamingEmitter emitter = new StreamingEmitter(lib);
            emitter.Devices.Add(device);
            emitter.SetControlPointCount(1, AdjustRate.None);

            int sc = 0; DateTimeOffset t0 = DateTimeOffset.MinValue, t1 = DateTimeOffset.MinValue;
            emitter.EmissionCallback = (em, iv, dl) => { foreach (var s in iv) { if (t0 == DateTimeOffset.MinValue) t0 = s.Time; t1 = s.Time; s.Points[0].Position = new Vector3(0f, 0f, 0.17f); s.Points[0].Intensity = 0f; sc++; } };
            emitter.Start(); Thread.Sleep(1000); emitter.Stop(); emitter.EmissionCallback = null;
            Console.WriteLine($"OK ({sc / (t1 - t0).TotalSeconds:F0} Hz)\n");

            Console.Write("  Participant ID: ");
            int pid; if (!int.TryParse(Console.ReadLine(), out pid)) pid = 999;

            try { new PairedComparisonExperiment(emitter).Run(pid); }
            catch (Exception ex) { Console.WriteLine($"\n  ERROR: {ex.Message}\n{ex.StackTrace}"); }
            finally { emitter.Stop(); emitter.Devices.Remove(device); }
        }

        static void PrintHeader(string title) => Console.WriteLine($"\n{title}\n{new string('-', title.Length)}\n");
    }
}