
import java.awt.event.ActionEvent;
import java.net.URL;
import java.util.ResourceBundle;
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.XYChart.Series;
import javafx.scene.chart.XYChart.Data;
import javafx.scene.control.TextField;
import javafx.scene.control.CheckBox;

//2yx-y^2+5-x^2   x0 = 0    y0 = 1   X = 3
public class FXMLDocumentController implements Initializable {

    @FXML
    private LineChart<Number, Number> MyChart;
    @FXML
    private LineChart<Number, Number> ErrorChart;
    @FXML
    private LineChart<Number, Number> GridChart;
    @FXML
    private TextField x0value;
    @FXML
    private TextField y0value;
    @FXML
    private TextField Xvalue;
    @FXML
    private TextField Nvalue;
    @FXML
    private TextField n0value;
    @FXML
    private TextField NvalueGrid;

    @FXML
    private CheckBox exact;
    @FXML
    private CheckBox euler;
    @FXML
    private CheckBox impEuler;
    @FXML
    private CheckBox rungeKutta;
    @FXML
    private CheckBox eulerGrid;
    @FXML
    private CheckBox impEulerGrid;
    @FXML
    private CheckBox rungeKuttaGrid;

    @FXML
    private void handleAction() {
        if (!x0value.getText().equals("")){
            x0 = Double.parseDouble(x0value.getText());
        }
        if (!y0value.getText().equals("")){
            y0 = Double.parseDouble(y0value.getText());
        }
        if (!Xvalue.getText().equals("")){
            X = Double.parseDouble(Xvalue.getText());
        }
        if (!Nvalue.getText().equals("")){
            N = Integer.parseInt(Nvalue.getText());
        }
        MyChart.getData().removeAll(Exact, Euler, ImpEuler, RungeKutta);
        ErrorChart.getData().removeAll(EulerError, ImprovedEulerError, RungeKuttaError);
        if (exact.isSelected()){
            ExactGraph(x0,y0,X,N);
        }
        if (euler.isSelected()){
            EulerGraph(x0,y0,X,N);
            EulerError(N);
        }
        if (impEuler.isSelected()){
            ImprovedEulerGraph(x0,y0,X,N);
            ImprovedEulerError(N);
        }
        if (rungeKutta.isSelected()){
            RungeKuttaGraph(x0,y0,X,N);
            RungeKuttaError(N);
        }
    }

    @FXML
    public void handleActionGrid(){
        if (!n0value.getText().equals("")){
            n0 = Integer.parseInt(n0value.getText());
        }
        if (!NvalueGrid.getText().equals("")){
            N = Integer.parseInt(NvalueGrid.getText());
        }
        GridChart.getData().removeAll(GridEuler, GridImprovedEuler, GridRungeKutta);
        if(eulerGrid.isSelected()){
            gridErrorEuler(n0, N);
        }
        if(impEulerGrid.isSelected()){
            gridErrorImprovedEuler(n0, N);
        }
        if(rungeKuttaGrid.isSelected()){
            gridErrorRungeKutta(n0, N);
        }
    }

    private Series Euler;
    private Series ImpEuler;
    private Series RungeKutta;
    private Series Exact;
    private Series EulerError;
    private Series ImprovedEulerError;
    private Series RungeKuttaError;
    private Series GridEuler;
    private Series GridImprovedEuler;
    private Series GridRungeKutta;
    private double e_x[], e_y[], ie_x[], ie_y[], rk_x[], rk_y[], ex_x[], ex_y[];
    private double x0 = 0, y0 = 1, X = 3;
    private int n0 = 1, N = 50;


    public double f(double x, double y){
        return 2*y*x - y*y + 5 - x*x;
    }

    @Override
    public void initialize(URL url, ResourceBundle rb) {

        EulerGraph(x0,y0,X,N);
        ImprovedEulerGraph(x0,y0,X,N);
        RungeKuttaGraph(x0,y0,X,N);
        ExactGraph(x0,y0,X,N);
        EulerError(N);
        ImprovedEulerError(N);
        RungeKuttaError(N);
        gridErrorEuler(n0, N);
        gridErrorImprovedEuler(n0, N);
        gridErrorRungeKutta(n0, N);
    }

    public void CountEulerGraph(double x0, double y0, double X, int N){
        e_x = new double[N];
        e_y = new double[N];
        double h = (X - x0) / N;

        e_x[0]=x0;
        e_y[0]=y0;
        for (int i = 1; i < N; i++){
            e_x[i] += i*h;
            e_y[i] = e_y[i-1] + h*f(e_x[i-1], e_y[i-1]);
        }
    }

    public void EulerGraph(double x0, double y0, double X, int N){
        CountEulerGraph(x0, y0, X, N);
        Euler = new Series();
        Euler.setName("Euler");
        for (int i=0; i<e_x.length; i++){
            Euler.getData().add(new Data(e_x[i], e_y[i]));
        }
        MyChart.getData().addAll(Euler);
    }

    public void CountImprovedEulerGraph(double x0, double y0, double X, int N){
        ie_x = new double[N];
        ie_y = new double[N];
        double h = (X - x0) / N;
        ie_x[0]=x0;
        ie_y[0]=y0;
        for (int i = 1; i < N; i++){
            ie_x[i] += i*h;
            ie_y[i] = ie_y[i-1] + (h/2)*(f(ie_x[i-1], ie_y[i-1]) +
                    f(ie_x[i], (ie_y[i-1] + h*f(ie_x[i-1], ie_y[i-1])))) ;
        }
    }

    public void ImprovedEulerGraph(double x0, double y0, double X, int N){
        CountImprovedEulerGraph(x0, y0, X, N);
        ImpEuler = new Series();
        ImpEuler.setName("Improved Euler");
        for (int i=0; i<ie_x.length; i++){
            ImpEuler.getData().add(new Data(ie_x[i], ie_y[i]));
        }
        MyChart.getData().addAll(ImpEuler);
    }

    public void CountRungeKuttaGraph(double x0, double y0, double X, int N){
        rk_x = new double[N];
        rk_y = new double[N];
        double h = (X - x0) / N;
        rk_x[0]=x0;
        rk_y[0]=y0;
        for (int i = 1; i < N; i++){
            rk_x[i] += i*h;
            rk_y[i] = rk_y[i-1] + (h/6)*(f(rk_x[i-1], rk_y[i-1]) +
                    2*f(rk_x[i-1] + h/2, rk_y[i-1] + (h/2)*f(rk_x[i-1], rk_y[i-1])) +
                    2*f(rk_x[i-1] + h/2, rk_y[i-1] + (h/2)*f(rk_x[i-1] + h/2, rk_y[i-1] +
                            (h/2)*f(rk_x[i-1], rk_y[i-1]))) +
                    f(rk_x[i-1] + h, rk_y[i-1] +
                            h*f(rk_x[i-1] + h/2, rk_y[i-1] + (h/2)*f(rk_x[i-1] + h/2, rk_y[i-1] +
                                    (h/2)*f(rk_x[i-1], rk_y[i-1])))));
        }
    }

    public void RungeKuttaGraph(double x0, double y0, double X, int N){
        CountRungeKuttaGraph(x0, y0, X, N);
        RungeKutta = new Series();
        RungeKutta.setName("Runge-Kutta");
        for (int i=0; i<rk_x.length; i++){
            RungeKutta.getData().add(new Data(rk_x[i], rk_y[i]));
        }
        MyChart.getData().addAll(RungeKutta);
    }

    public void CountExactGraph(double x0, double y0, double X, int N){
        ex_x = new double[N];
        ex_y = new double[N];
        double h = (X - x0) / N;
        ex_x[0]=x0;
        ex_y[0]=y0;
        double C = (1/(y0 - x0 - 2.0) + 1.0/4.0)/Math.exp(4*x0);
        for (int i = 1; i < N; i++){
            ex_x[i] += i*h;
            ex_y[i] = 1/(C*Math.exp(4*ex_x[i]) - 1.0/4.0) + ex_x[i] + 2.0;
        }
    }

    public void ExactGraph(double x0, double y0, double X, int N){
        CountExactGraph(x0, y0, X, N);
        Exact = new Series();
        Exact.setName("Exact");
        for (int i = 0; i<ex_x.length; i++){
            Exact.getData().add(new Data(ex_x[i], ex_y[i]));
        }
        MyChart.getData().addAll(Exact);
    }

    public void EulerError(int N){
        double er_x[] = new double[N];
        double er_y[] = new double[N];
        for (int i = 0; i < e_x.length; i++){
            er_x[i] = e_x[i];
            er_y[i] = Math.abs(e_y[i] - ex_y[i]);
        }
        EulerError = new Series();
        EulerError.setName("Euler Error");
        for (int i = 0; i<er_x.length; i++){
            EulerError.getData().add(new Data(er_x[i], er_y[i]));
        }
        ErrorChart.getData().addAll(EulerError);
    }

    public void ImprovedEulerError(int N){
        double i_er_x[] = new double[N];
        double i_er_y[] = new double[N];
        for (int i = 0; i < ie_x.length; i++){
            i_er_x[i] = ie_x[i];
            i_er_y[i] = Math.abs(ie_y[i] - ex_y[i]);
        }
        ImprovedEulerError = new Series();
        ImprovedEulerError.setName("Imp.Euler Error");
        for (int i = 0; i<i_er_x.length; i++){
            ImprovedEulerError.getData().add(new Data(i_er_x[i], i_er_y[i]));
        }
        ErrorChart.getData().addAll(ImprovedEulerError);
    }

    public void RungeKuttaError(int N){
        double er_x[] = new double[N];
        double er_y[] = new double[N];
        for (int i = 0; i < rk_x.length; i++){
            er_x[i] = rk_x[i];
            er_y[i] = Math.abs(rk_y[i] - ex_y[i]);
        }
        RungeKuttaError = new Series();
        RungeKuttaError.setName("Runge-Kutta Error");
        for (int i = 0; i<er_x.length; i++){
            RungeKuttaError.getData().add(new Data(er_x[i], er_y[i]));
        }
        ErrorChart.getData().addAll(RungeKuttaError);
    }

    public void gridErrorEuler(int n0, int N){
        double gr_e_x[] = new double[N-n0];
        double gr_e_y[] = new double[N-n0];
        double max = 0;
        for (int j = n0; j<N; j++){
            CountEulerGraph(x0, y0, X, j);
            CountExactGraph(x0, y0, X, j);
            for (int k = 0; k<e_x.length; k++) {
                if (Math.abs(e_y[k] - ex_y[k]) > max) {
                    max = Math.abs(e_y[k] - ex_y[k]);
                }
            }
            gr_e_x[j-n0] = j;
            gr_e_y[j-n0] = max;
            max = 0;
        }
        GridEuler = new Series();
        GridEuler.setName("Grid Euler");
        for (int i = 0; i<gr_e_x.length; i++){
            GridEuler.getData().add(new Data(gr_e_x[i], gr_e_y[i]));
        }
        GridChart.getData().addAll(GridEuler);
    }

    public void gridErrorImprovedEuler(int n0, int N){
        double gr_ie_x[] = new double[N-n0];
        double gr_ie_y[] = new double[N-n0];
        double max = 0;
        for (int j = n0; j<N; j++){
            CountImprovedEulerGraph(x0, y0, X, j);
            CountExactGraph(x0, y0, X, j);
            for (int k = 0; k<ie_x.length; k++) {
                if (Math.abs(ie_y[k] - ex_y[k]) > max) {
                    max = Math.abs(ie_y[k] - ex_y[k]);
                }
            }
            gr_ie_x[j-n0] = j;
            gr_ie_y[j-n0] = max;
            max = 0;
        }
        GridImprovedEuler = new Series();
        GridImprovedEuler.setName("Grid Improved Euler");
        for (int i = 0; i<gr_ie_x.length; i++){
            GridImprovedEuler.getData().add(new Data(gr_ie_x[i], gr_ie_y[i]));
        }
        GridChart.getData().addAll(GridImprovedEuler);
    }

    public void gridErrorRungeKutta(int n0, int N){
        double gr_rk_x[] = new double[N-n0];
        double gr_rk_y[] = new double[N-n0];
        double max = 0;
        for (int j = n0; j<N; j++){
            CountRungeKuttaGraph(x0, y0, X, j);
            CountExactGraph(x0, y0, X, j);
            for (int k = 0; k<rk_x.length; k++) {
                if (Math.abs(rk_y[k] - ex_y[k]) > max) {
                    max = Math.abs(rk_y[k] - ex_y[k]);
                }
            }
            gr_rk_x[j-n0] = j;
            gr_rk_y[j-n0] = max;
            max = 0;
        }
        GridRungeKutta = new Series();
        GridRungeKutta.setName("Grid RungeKutta");
        for (int i = 0; i<gr_rk_x.length; i++){
            GridRungeKutta.getData().add(new Data(gr_rk_x[i], gr_rk_y[i]));
        }
        GridChart.getData().addAll(GridRungeKutta);
    }
}
